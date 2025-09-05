#!/usr/bin/env python3

import pysam
import sys
import os
import re
from collections import defaultdict, OrderedDict
import logging

def setup_logging():
    """设置日志"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

# --- 函数定义 (get_caller_priority_and_name, get_caller_display_name, get_input_files 保持不变) ---
def clean_chromosome_name(chrom_str):
    """清理染色体名称，处理VarScan2的特殊格式"""
    chrom_str = str(chrom_str).strip()
    
    # 处理VarScan2的特殊格式: "dna/variants/varscan2/PATIENT1.indel.vcf:chr6" -> "chr6"
    if ':chr' in chrom_str:
        chrom_part = chrom_str.split(':chr')[-1]
        return f"chr{chrom_part}"
    elif ':' in chrom_str and chrom_str.split(':')[-1].isdigit():
        # 处理 "path:1" -> "chr1" 的情况
        chrom_part = chrom_str.split(':')[-1]
        return f"chr{chrom_part}"
    elif ':' in chrom_str and chrom_str.split(':')[-1] in ['X', 'Y', 'M', 'MT']:
        # 处理 "path:X" -> "chrX" 的情况
        chrom_part = chrom_str.split(':')[-1]
        return f"chr{chrom_part}"
    elif chrom_str.startswith('chr'):
        return chrom_str
    else:
        # 如果是纯数字或X/Y，添加chr前缀
        if chrom_str.isdigit() or chrom_str in ['X', 'Y', 'M', 'MT']:
            return f"chr{chrom_str}"
        else:
            return chrom_str

def normalize_variant_key(record):
    """标准化变异位点键值，用于比较不同caller的结果"""
    # 从record的contig获取原始染色体名称
    original_chrom = str(record.contig)
    clean_chrom = clean_chromosome_name(original_chrom)
    
    # 使用清理后的染色体、位置、参考碱基、替代碱基作为键
    alts = []
    if record.alts:
        alts = [str(alt) for alt in record.alts if alt is not None]
    
    key = (
        clean_chrom,
        int(record.start),  # 0-based position
        str(record.ref),
        tuple(sorted(alts))
    )
    return key

def get_caller_priority_and_name(filepath):
    """从文件路径获取caller名称和优先级"""
    basename = os.path.basename(str(filepath))
    
    if 'mutect2' in basename.lower() or 'mutect' in basename.lower():
        return 1, 'mutect2'
    elif 'strelka' in basename.lower():
        return 2, 'strelka2'
    elif 'varscan' in basename.lower():
        return 3, 'varscan2'
    elif 'manta' in basename.lower():
        return 4, 'manta'
    else:
        return 999, 'unknown'

def get_caller_display_name(caller_name):
    """获取caller的显示名称"""
    name_map = {
        'mutect2': 'MuTect2',
        'strelka2': 'Strelka2', 
        'varscan2': 'VarScan2',
        'manta': 'Manta'
    }
    return name_map.get(caller_name, caller_name.capitalize())

def get_input_files(input_obj):
    """从Snakemake输入对象中提取文件列表"""
    files = []
    
    try:
        # 尝试直接迭代
        for item in input_obj:
            filepath = str(item)
            if filepath and os.path.exists(filepath):
                files.append(filepath)
        return files
    except:
        pass
    
    try:
        # 尝试作为具名输入处理
        if hasattr(input_obj, '_names'):
            for name in input_obj._names:
                filepath = str(getattr(input_obj, name))
                if filepath and os.path.exists(filepath):
                    files.append(filepath)
        elif hasattr(input_obj, '__dict__'):
            for key, value in input_obj.__dict__.items():
                if not key.startswith('_'):
                    filepath = str(value)
                    if filepath and os.path.exists(filepath):
                        files.append(filepath)
        return files
    except:
        pass
    
    # 最后尝试转为字符串
    try:
        filepath = str(input_obj)
        if filepath and os.path.exists(filepath):
            files.append(filepath)
    except:
        pass
    
    return files


def build_sample_map(sample_id):
    """# 核心修复: 创建样本名映射"""
    return {
        'NORMAL': f"{sample_id}_normal",
        'TUMOR': f"{sample_id}_tumor",
        f"{sample_id}_normal": f"{sample_id}_normal",
        f"{sample_id}_tumor": f"{sample_id}_tumor"
    }

def merge_headers(file_list, sample_map):
    """
    **修复版**: 合并所有VCF文件的header信息。
    通过创建一个全新的header对象来避免API问题，并正确处理样本名。
    """
    logger = logging.getLogger(__name__)
    
    # 收集所有header信息
    all_info_fields = {}
    all_format_fields = {}
    all_filter_fields = {}
    all_contigs = {} # 使用字典来存储contig和其长度
    all_samples = set() # 存储规范化后的样本名
    template_header = None # 仍然需要一个模板来获取基本信息
    
    for filepath in file_list:
        try:
            with pysam.VariantFile(filepath) as vcf:
                logger.info(f"合并header from: {filepath}")
                
                if template_header is None:
                    template_header = vcf.header.copy()
                
                for info_id, info_record in vcf.header.info.items():
                    all_info_fields[info_id] = info_record
                
                for format_id, format_record in vcf.header.formats.items():
                    all_format_fields[format_id] = format_record
                
                for filter_id, filter_record in vcf.header.filters.items():
                    all_filter_fields[filter_id] = filter_record
                
                for contig_id, contig_record in vcf.header.contigs.items():
                    all_contigs[contig_id] = contig_record
                
                # 使用样本名映射
                for sample in vcf.header.samples:
                    all_samples.add(sample_map.get(sample, sample))
                    
        except Exception as e:
            logger.error(f"读取header时出错 {filepath}: {str(e)}")
            continue
    
    if template_header is None:
        raise ValueError("没有有效的VCF文件")
    
    # **核心修复**: 创建一个全新的、空的header，而不是修改模板
    merged_header = pysam.VariantHeader()
    
    # 1. 添加所有contigs
    for contig_id, contig_record in all_contigs.items():
        try:
            merged_header.contigs.add(contig_id, length=contig_record.length)
        except Exception as e:
            logger.warning(f"无法添加contig {contig_id}: {e}")

    # 2. 添加所有INFO字段
    for info_id, info_record in all_info_fields.items():
        try:
            merged_header.info.add(info_id, info_record.number, info_record.type, info_record.description)
        except Exception as e:
            logger.warning(f"无法添加INFO字段 {info_id}: {e}")

    # 3. 添加所有FORMAT字段
    for format_id, format_record in all_format_fields.items():
        try:
            merged_header.formats.add(format_id, format_record.number, format_record.type, format_record.description)
        except Exception as e:
            logger.warning(f"无法添加FORMAT字段 {format_id}: {e}")

    # 4. 添加所有FILTER字段
    for filter_id, filter_record in all_filter_fields.items():
        try:
            merged_header.filters.add(filter_id, description=filter_record.description)
        except Exception as e:
            logger.warning(f"无法添加FILTER字段 {filter_id}: {e}")
            
    # 5. 添加自定义的合并字段
    try:
        merged_header.info.add('CALLERS', '.', 'String', 'Variant callers that detected this variant')
        merged_header.info.add('NCALLERS', 1, 'Integer', 'Number of callers that detected this variant')
        merged_header.info.add('CALLER_SUPPORT', '.', 'String', 'Caller support summary (detected/total: caller_list)')
    except Exception as e:
        logger.warning(f"无法添加合并字段: {e}")

    # 6. 最后，添加规范化、排序后的样本名
    for sample_name in sorted(list(all_samples)):
        merged_header.samples.add(sample_name)
    
    logger.info(f"Header合并完成:")
    logger.info(f"  样本数: {len(merged_header.samples)} - {list(merged_header.samples)}")
    
    return merged_header

def merge_variant_records(caller_records, merged_header, sample_map):
    """合并同一位点的多个变异记录"""
    logger = logging.getLogger(__name__)
    
    caller_records.sort(key=lambda x: x[0])
    
    primary_priority, primary_caller, primary_record, clean_chrom = caller_records[0]
    
    variant_info = {
        'chrom': clean_chrom, 'pos': primary_record.start, 'stop': primary_record.stop,
        'id': primary_record.id, 'ref': primary_record.ref, 'alts': primary_record.alts,
        'qual': primary_record.qual, 'filters': set(), 'info': {}, 'samples': {}
    }
    
    detected_callers = []
    
    for priority, caller_name, record, chrom in caller_records:
        if caller_name not in detected_callers:
            detected_callers.append(caller_name)
        
        for filter_name in record.filter:
            variant_info['filters'].add(filter_name)
        
        for key, value in record.info.items():
            if key not in ['CALLERS', 'NCALLERS', 'CALLER_SUPPORT']:
                if key not in variant_info['info'] or caller_name == primary_caller:
                     variant_info['info'][key] = value
        
        # # 核心修复: 使用样本名映射
        for sample_name in record.samples:
            norm_sample_name = sample_map.get(sample_name, sample_name)
            if norm_sample_name not in variant_info['samples']:
                variant_info['samples'][norm_sample_name] = {}
            
            sample_data = record.samples[sample_name]
            for format_key, format_value in sample_data.items():
                if format_key not in variant_info['samples'][norm_sample_name] or caller_name == primary_caller:
                    variant_info['samples'][norm_sample_name][format_key] = format_value
    
    return variant_info, detected_callers

def create_output_record(variant_info, detected_callers, all_available_callers, vcf_out):
    """创建输出VCF记录"""
    logger = logging.getLogger(__name__)
    
    try:
        new_record = vcf_out.new_record(
            contig=variant_info['chrom'], start=variant_info['pos'], stop=variant_info['stop'],
            alleles=[variant_info['ref']] + list(variant_info['alts']) if variant_info['alts'] else [variant_info['ref']]
        )
        
        if variant_info['id']: new_record.id = variant_info['id']
        if variant_info['qual'] is not None: new_record.qual = variant_info['qual']
        
        if variant_info['filters']:
            for filter_name in variant_info['filters']:
                if filter_name != 'PASS': new_record.filter.add(filter_name)
        else:
            new_record.filter.add('PASS')
        
        # # 核心修复: 智能处理多等位基因INFO字段
        num_alts = len(new_record.alts)
        for key, value in variant_info['info'].items():
            try:
                info_field_meta = vcf_out.header.info.get(key)
                if info_field_meta:
                    number = info_field_meta.number
                    # 如果字段是 per-alternate-allele ('A') 或 per-allele ('R')
                    if number in ('A', 'R') and isinstance(value, (list, tuple)):
                        # 确保值的数量与ALT数量匹配
                        if len(value) != num_alts:
                            logger.debug(f"调整INFO字段 {key} 的值数量以匹配 {num_alts} 个ALT。原始值: {value}")
                            # 简单的填充/截断策略
                            value = tuple(list(value)[:num_alts] + ['.'] * (num_alts - len(value)))
                new_record.info[key] = value
            except Exception as e:
                logger.warning(f"设置INFO字段 {key} (值: {value}) 失败: {e}")
        
        try:
            detected_callers_clean = [c.strip() for c in detected_callers]
            caller_display_names = [get_caller_display_name(name) for name in detected_callers_clean]
            new_record.info['CALLERS'] = ','.join(caller_display_names)
            new_record.info['NCALLERS'] = len(detected_callers)
            
            total_possible = len(all_available_callers)
            support_info = f"{len(detected_callers)}/{total_possible}:{','.join(caller_display_names)}"
            new_record.info['CALLER_SUPPORT'] = support_info
        except Exception as e:
            logger.warning(f"添加caller信息失败: {e}")
        
        for header_sample in vcf_out.header.samples:
            if header_sample in variant_info['samples']:
                sample_data = variant_info['samples'][header_sample]
                for format_key, format_value in sample_data.items():
                    try:
                        new_record.samples[header_sample][format_key] = format_value
                    except Exception as e:
                        logger.debug(f"设置样本 {header_sample} 的字段 {format_key} 失败: {e}")
            else:
                if 'GT' in vcf_out.header.formats:
                    new_record.samples[header_sample]['GT'] = (None, None)
                        
        return new_record
        
    except Exception as e:
        logger.error(f"创建输出记录失败: {e}")
        import traceback
        logger.error(f"详细错误: {traceback.format_exc()}")
        return None

def merge_variants(input_files, output_file, sample_id):
    """合并多个caller的变异结果"""
    logger = setup_logging()
    logger.info(f"开始合并变异文件: {input_files}")
    logger.info(f"输出文件: {output_file}")
    
    file_list = get_input_files(input_files)
    if not file_list:
        logger.error("没有找到有效的输入文件"); sys.exit(1)
    
    # # 核心修复: 创建并使用样本名映射
    sample_map = build_sample_map(sample_id)
    
    prioritized_files = []
    all_available_callers = []
    for filepath in file_list:
        priority, caller_name = get_caller_priority_and_name(filepath)
        prioritized_files.append((priority, caller_name, filepath))
        all_available_callers.append(caller_name)
    prioritized_files.sort()
    
    logger.info("开始合并VCF headers...")
    merged_header = merge_headers(file_list, sample_map)
    
    variants_by_position = defaultdict(list)
    
    for priority, caller_name, filepath in prioritized_files:
        display_name = get_caller_display_name(caller_name)
        logger.info(f"处理 {display_name} 的结果: {filepath}")
        try:
            with pysam.VariantFile(filepath) as vcf_in:
                for record in vcf_in:
                    pos_key = normalize_variant_key(record)
                    variants_by_position[pos_key].append((priority, caller_name, record, pos_key[0]))
        except Exception as e:
            logger.error(f"读取文件 {filepath} 时出错: {e}")
            continue
    
    logger.info(f"总共发现 {len(variants_by_position)} 个唯一变异位点")
    
    with pysam.VariantFile(output_file, 'wz', header=merged_header) as vcf_out:
        sorted_variants = sorted(variants_by_position.items(), key=lambda x: (
            (0, int(x[0][0][3:])) if x[0][0].startswith('chr') and x[0][0][3:].isdigit()
            else (1, 0) if x[0][0] == 'chrX' else (2, 0) if x[0][0] == 'chrY'
            else (3, 0) if x[0][0] in ['chrM', 'chrMT'] else (4, x[0][0]),
            x[0][1]
        ))
        
        for pos_key, caller_records in sorted_variants:
            try:
                variant_info, detected_callers = merge_variant_records(caller_records, merged_header, sample_map)
                new_record = create_output_record(variant_info, detected_callers, all_available_callers, vcf_out)
                if new_record: vcf_out.write(new_record)
            except Exception as e:
                logger.error(f"处理变异 {pos_key} 时出错: {e}")
                continue

    try:
        logger.info("创建tabix索引...")
        pysam.tabix_index(output_file, preset='vcf', force=True)
        logger.info("索引创建完成")
    except Exception as e:
        logger.error(f"创建索引时出错: {e}")

    logger.info("合并完成")

def main():
    """主函数 - 命令行模式"""
    if len(sys.argv) < 4:
        print("用法: python merge_variants.py <sample_id> <output_file> <input1.vcf> [<input2.vcf> ...]")
        sys.exit(1)
    
    sample_id = sys.argv[1]
    output_file = sys.argv[2]
    input_files = sys.argv[3:]
    
    merge_variants(input_files, output_file, sample_id)

if __name__ == "__main__":
    if 'snakemake' in globals():
        merge_variants(snakemake.input, snakemake.output.vcf, snakemake.params.sample_id)
    else:
        main()