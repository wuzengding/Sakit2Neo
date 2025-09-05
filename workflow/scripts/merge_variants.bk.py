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
    
    #  键只包含染色体和起始位置
    key = (
        clean_chrom,
        int(record.start),  # 0-based position
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

def merge_headers(file_list):
    """合并所有VCF文件的header信息"""
    logger = logging.getLogger(__name__)
    
    # 收集所有header信息
    all_info_fields = {}
    all_format_fields = {}
    all_filter_fields = {}
    all_contigs = set()
    all_samples = set()
    template_header = None
    
    for filepath in file_list:
        try:
            with pysam.VariantFile(filepath) as vcf:
                logger.info(f"合并header from: {filepath}")
                
                # 保存第一个文件的header作为模板
                if template_header is None:
                    template_header = vcf.header.copy()
                
                # 收集INFO字段
                for info_id in vcf.header.info:
                    info_record = vcf.header.info[info_id]
                    all_info_fields[info_id] = info_record
                
                # 收集FORMAT字段
                for format_id in vcf.header.formats:
                    format_record = vcf.header.formats[format_id]
                    all_format_fields[format_id] = format_record
                
                # 收集FILTER字段
                for filter_id in vcf.header.filters:
                    filter_record = vcf.header.filters[filter_id]
                    all_filter_fields[filter_id] = filter_record
                
                # 收集contigs
                for contig_id in vcf.header.contigs:
                    all_contigs.add(contig_id)
                
                # 收集样本名
                for sample in vcf.header.samples:
                    all_samples.add(sample)
                    
        except Exception as e:
            logger.error(f"读取header时出错 {filepath}: {str(e)}")
            continue
    
    if template_header is None:
        raise ValueError("没有有效的VCF文件")
    
    # 创建合并后的header
    merged_header = template_header.copy()
    
    # 添加所有INFO字段
    for info_id, info_record in all_info_fields.items():
        if info_id not in merged_header.info:
            try:
                merged_header.info.add(
                    info_id,
                    number=info_record.number,
                    type=info_record.type,
                    description=info_record.description
                )
                logger.debug(f"添加INFO字段: {info_id}")
            except Exception as e:
                logger.warning(f"无法添加INFO字段 {info_id}: {e}")
    
    # 添加所有FORMAT字段
    for format_id, format_record in all_format_fields.items():
        if format_id not in merged_header.formats:
            try:
                merged_header.formats.add(
                    format_id,
                    number=format_record.number,
                    type=format_record.type,
                    description=format_record.description
                )
                logger.debug(f"添加FORMAT字段: {format_id}")
            except Exception as e:
                logger.warning(f"无法添加FORMAT字段 {format_id}: {e}")
    
    # 添加所有FILTER字段
    for filter_id, filter_record in all_filter_fields.items():
        if filter_id not in merged_header.filters:
            try:
                # 修复FILTER字段添加方法
                merged_header.filters.add(filter_id, description=filter_record.description)
                logger.debug(f"添加FILTER字段: {filter_id}")
            except Exception as e:
                logger.warning(f"无法添加FILTER字段 {filter_id}: {e}")
    
    # 添加所有contigs
    for contig_id in all_contigs:
        if contig_id not in merged_header.contigs:
            try:
                merged_header.contigs.add(contig_id)
                logger.debug(f"添加contig: {contig_id}")
            except Exception as e:
                logger.warning(f"无法添加contig {contig_id}: {e}")
    
    # 添加合并特有的字段
    try:
        if 'CALLERS' not in merged_header.info:
            merged_header.info.add('CALLERS', number='.', type='String', 
                                 description='Variant callers that detected this variant')
        if 'NCALLERS' not in merged_header.info:
            merged_header.info.add('NCALLERS', number=1, type='Integer', 
                                 description='Number of callers that detected this variant')
        if 'CALLER_SUPPORT' not in merged_header.info:
            merged_header.info.add('CALLER_SUPPORT', number='.', type='String', 
                                 description='Caller support summary (detected/total: caller_list)')
    except Exception as e:
        logger.warning(f"无法添加合并字段: {e}")
    
    # 确保所有样本都在header中
    for sample_name in all_samples:
        if sample_name not in merged_header.samples:
            try:
                merged_header.samples.add(sample_name)
                logger.debug(f"添加样本到header: {sample_name}")
            except Exception as e:
                logger.warning(f"无法添加样本 {sample_name}: {e}")
    
    logger.info(f"Header合并完成:")
    logger.info(f"  INFO字段数: {len(merged_header.info)}")
    logger.info(f"  FORMAT字段数: {len(merged_header.formats)}")
    logger.info(f"  FILTER字段数: {len(merged_header.filters)}")
    logger.info(f"  Contigs数: {len(merged_header.contigs)}")
    logger.info(f"  样本数: {len(all_samples)} - {list(all_samples)}")
    
    return merged_header, all_samples


def merge_variant_records(caller_records, merged_header):
    """
    **修改版**: 合并同一【位置】的多个变异记录。
    它会合并所有独特的等位基因，并以优先级最高的记录为基础来填充元信息。
    """
    logger = logging.getLogger(__name__)
    
    # 按优先级排序（mutect2优先）
    caller_records.sort(key=lambda x: x[0])
    
    # 选择优先级最高的记录作为主记录，它的元信息（如ID, QUAL, FILTERs, INFO, SAMPLEs）将被优先使用
    primary_priority, primary_caller, primary_record, clean_chrom = caller_records[0]
    
    # --- 新增核心逻辑: 收集并合并所有等位基因 ---
    all_refs = set()
    all_alts = set()
    detected_callers = []
    
    for priority, caller_name, record, chrom in caller_records:
        all_refs.add(record.ref)
        if record.alts:
            for alt in record.alts:
                all_alts.add(alt)
        
        # 收集独特的caller名称
        if caller_name not in detected_callers:
            detected_callers.append(caller_name)

    # 确定最终的REF和ALT
    # VCF规范要求一个位置只有一个REF。如果存在多个REF（例如，相邻INDEL规范化后），
    # 我们选择其中最短的一个作为公共REF，并相应调整其他等位基因。
    # 这是一个健壮的处理方式。
    
    if len(all_refs) > 1:
        logger.warning(f"多个REF在同一位置 {clean_chrom}:{primary_record.start + 1}. "
                       f"Refs: {all_refs}. 将尝试寻找公共REF。")
        
        # 找到最短的REF作为公共REF
        final_ref = min(all_refs, key=len)
        
        # 调整所有的REF和ALT
        adjusted_alts = set()
        for ref in all_refs:
            if ref.startswith(final_ref):
                suffix = ref[len(final_ref):]
                if not suffix: # 这是final_ref本身
                    adjusted_alts.add('.') # '.'表示REF和ALT相同，后面会过滤掉
                else:
                    adjusted_alts.add(suffix) # 例如 REF=AG, final_ref=A -> 添加 G
            else:
                 logger.error(f"无法处理不兼容的REF: {ref} vs {final_ref}")
                 # 在这种复杂情况下，可以简单地跳过，或者采取其他策略
                 # 这里我们简单地将它们都加入，让下游工具处理
                 adjusted_alts.add(ref)

        for alt in all_alts:
            if alt.startswith(final_ref):
                suffix = alt[len(final_ref):]
                if not suffix:
                     adjusted_alts.add('.')
                else:
                    adjusted_alts.add(suffix)
            else:
                logger.error(f"无法处理不兼容的ALT: {alt} vs {final_ref}")
                adjusted_alts.add(alt)

        all_alts = adjusted_alts

    else:
        final_ref = list(all_refs)[0]

    # 从合并后的ALT中移除代表REF的'.'和任何空字符串
    final_alts = sorted([alt for alt in all_alts if alt and alt != '.'])
    
    # 如果所有ALT都被移除后列表为空，说明所有变异都是REF=ALT，这不应该发生
    # 但作为安全检查，我们可以保留主记录的ALT
    if not final_alts and primary_record.alts:
        final_alts = list(primary_record.alts)
    # --- 结束新增核心逻辑 ---


    # 创建新记录的基本信息，使用合并后的等位基因
    variant_info = {
        'chrom': clean_chrom,
        'pos': primary_record.start,
        'stop': primary_record.stop,
        'id': primary_record.id,
        'ref': final_ref,
        'alts': tuple(final_alts),
        'qual': primary_record.qual,
        'filters': set(),
        'info': {},
        'samples': {}
    }
    
    # 收集所有caller的元信息（FILTERs, INFO, SAMPLEs）
    # 这部分逻辑与您原来的版本保持一致，以主caller的记录为准
    all_samples = set()
    
    for priority, caller_name, record, chrom in caller_records:
        # 合并FILTER字段
        for filter_name in record.filter:
            variant_info['filters'].add(filter_name)
        
        # 合并INFO字段 - 以主caller为准
        for key, value in record.info.items():
            if key not in ['CALLERS', 'NCALLERS', 'CALLER_SUPPORT']:
                if key not in variant_info['info'] or caller_name == primary_caller:
                    variant_info['info'][key] = value
        
        # 合并样本信息 - 以主caller为准
        for sample_name in record.samples:
            all_samples.add(sample_name)
            if sample_name not in variant_info['samples']:
                variant_info['samples'][sample_name] = {}
            
            sample_data = record.samples[sample_name]
            for format_key, format_value in sample_data.items():
                if format_key not in variant_info['samples'][sample_name] or caller_name == primary_caller:
                    variant_info['samples'][sample_name][format_key] = format_value
    
    return variant_info, detected_callers, all_samples
    
def create_output_record(variant_info, detected_callers, all_available_callers, vcf_out):
    """创建输出VCF记录"""
    logger = logging.getLogger(__name__)
    
    try:
        # 创建新的变异记录
        new_record = vcf_out.new_record(
            contig=variant_info['chrom'],
            start=variant_info['pos'],
            stop=variant_info['stop'],
            alleles=[variant_info['ref']] + list(variant_info['alts']) if variant_info['alts'] else [variant_info['ref']]
        )
        
        # 设置基本字段
        if variant_info['id']:
            new_record.id = variant_info['id']
        if variant_info['qual'] is not None:
            new_record.qual = variant_info['qual']
        
        # 设置FILTER字段
        if variant_info['filters']:
            for filter_name in variant_info['filters']:
                if filter_name != 'PASS':  # PASS是默认值，需要明确设置
                    new_record.filter.add(filter_name)
        else:
            # 如果没有filter，设置为PASS
            new_record.filter.add('PASS')
        
        # 设置INFO字段
        for key, value in variant_info['info'].items():
            try:
                new_record.info[key] = value
            except Exception as e:
                logger.warning(f"设置INFO字段 {key} 失败: {e}")
        
        # 添加caller信息
        try:
            caller_display_names = [get_caller_display_name(name) for name in detected_callers]
            new_record.info['CALLERS'] = ','.join(caller_display_names)
            new_record.info['NCALLERS'] = len(detected_callers)
            
            # 创建支持信息
            total_possible = len(all_available_callers)
            support_info = f"{len(detected_callers)}/{total_possible}:{','.join(caller_display_names)}"
            new_record.info['CALLER_SUPPORT'] = support_info
        except Exception as e:
            logger.warning(f"添加caller信息失败: {e}")
        
        # 设置样本FORMAT字段 - 修复样本名称匹配问题
        header_samples = list(vcf_out.header.samples)  # 获取header中的样本名称
        logger.debug(f"Header中的样本: {header_samples}")
        logger.debug(f"数据中的样本: {list(variant_info['samples'].keys())}")
        
        # 为所有header中的样本设置数据
        for header_sample in header_samples:
            try:
                # 查找匹配的样本数据
                sample_data_found = None
                
                # 首先尝试精确匹配
                if header_sample in variant_info['samples']:
                    sample_data_found = variant_info['samples'][header_sample]
                else:
                    # 尝试模糊匹配
                    for data_sample_name, data_sample_data in variant_info['samples'].items():
                        # 检查是否为同一个样本的不同命名
                        if (header_sample.lower() in data_sample_name.lower() or 
                            data_sample_name.lower() in header_sample.lower()):
                            sample_data_found = data_sample_data
                            logger.debug(f"样本名称匹配: {header_sample} -> {data_sample_name}")
                            break
                
                if sample_data_found:
                    # 设置找到的样本数据
                    for format_key, format_value in sample_data_found.items():
                        try:
                            new_record.samples[header_sample][format_key] = format_value
                        except Exception as e:
                            logger.debug(f"设置样本 {header_sample} 的字段 {format_key} 失败: {e}")
                else:
                    # 如果找不到匹配的样本数据，设置默认值
                    logger.debug(f"样本 {header_sample} 没有找到匹配数据，使用默认值")
                    # 设置基本的GT字段
                    if 'GT' in vcf_out.header.formats:
                        new_record.samples[header_sample]['GT'] = (None, None)
                        
            except Exception as e:
                logger.warning(f"设置样本 {header_sample} 数据失败: {e}")
                continue
        
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
    
    # 获取输入文件列表
    file_list = get_input_files(input_files)
    logger.info(f"解析到的文件列表: {file_list}")
    
    if not file_list:
        logger.error("没有找到有效的输入文件")
        sys.exit(1)
    
    # 按优先级排序文件
    prioritized_files = []
    all_available_callers = []
    for filepath in file_list:
        priority, caller_name = get_caller_priority_and_name(filepath)
        prioritized_files.append((priority, caller_name, filepath))
        all_available_callers.append(caller_name)
    
    prioritized_files.sort(key=lambda x: x[0])  # 按优先级排序
    logger.info(f"按优先级排序的文件: {[(p, c, f) for p, c, f in prioritized_files]}")
    logger.info(f"可用的callers: {all_available_callers}")
    
    # 合并所有VCF的header
    logger.info("开始合并VCF headers...")
    merged_header, all_samples_in_files = merge_headers(file_list)
    
    # 存储所有变异
    variants_by_position = defaultdict(list)  # position_key -> [(priority, caller_name, record, clean_chrom)]
    
    # 读取每个caller的结果
    for priority, caller_name, filepath in prioritized_files:
        display_name = get_caller_display_name(caller_name)
        logger.info(f"处理 {display_name} 的结果: {filepath}")
        
        try:
            with pysam.VariantFile(filepath) as vcf_in:
                variant_count = 0
                for record in vcf_in.fetch():
                    pos_key = normalize_variant_key(record)
                    clean_chrom = pos_key[0]
                    
                    variants_by_position[pos_key].append((priority, caller_name, record, clean_chrom))
                    variant_count += 1
                
                logger.info(f"{display_name} 包含 {variant_count} 个变异")
        
        except Exception as e:
            logger.error(f"读取文件 {filepath} 时出错: {str(e)}")
            continue
    
    logger.info(f"总共发现 {len(variants_by_position)} 个唯一变异位点")
    
    # 创建输出VCF
    try:
        # 使用压缩格式输出
        with pysam.VariantFile(output_file, 'wz', header=merged_header) as vcf_out:
            # 统计信息
            stats = {
                'total_variants': 0,
                'by_caller_count': defaultdict(int),
                'by_primary_caller': defaultdict(int),
                'skipped_variants': 0
            }
            
            # 按染色体和位置排序处理变异
            sorted_variants = sorted(
                variants_by_position.items(),
                key=lambda x: (
                    # 自定义染色体排序：数字染色体 -> X -> Y -> M -> 其他
                    (0, int(x[0][0][3:])) if x[0][0].startswith('chr') and x[0][0][3:].isdigit()
                    else (1, 0) if x[0][0] == 'chrX'
                    else (2, 0) if x[0][0] == 'chrY' 
                    else (3, 0) if x[0][0] in ['chrM', 'chrMT']
                    else (4, x[0][0]),
                    x[0][1]  # 按位置排序
                )
            )
            
            for pos_key, caller_records in sorted_variants:
                try:
                    # 合并同一位点的多个变异记录
                    variant_info, detected_callers, samples_for_variant = merge_variant_records(
                        caller_records, merged_header
                    )
                    
                    # 创建输出记录
                    new_record = create_output_record(
                        variant_info, detected_callers, all_available_callers, vcf_out
                    )
                    
                    if new_record is not None:
                        # 写入合并后的记录
                        vcf_out.write(new_record)
                        
                        # 更新统计信息
                        stats['total_variants'] += 1
                        caller_count = len(caller_records)
                        stats['by_caller_count'][caller_count] += 1
                        
                        primary_caller = caller_records[0][1]  # 已按优先级排序
                        stats['by_primary_caller'][primary_caller] += 1
                        
                        logger.debug(f"成功写入变异: {variant_info['chrom']}:{pos_key[1]+1} "
                                   f"({len(detected_callers)} callers: {','.join(detected_callers)})")
                    else:
                        stats['skipped_variants'] += 1
                        logger.warning(f"跳过变异: {pos_key}")
                
                except Exception as e:
                    logger.error(f"处理变异 {pos_key} 时出错: {str(e)}")
                    import traceback
                    logger.error(f"详细错误: {traceback.format_exc()}")
                    stats['skipped_variants'] += 1
                    continue
    
    except Exception as e:
        logger.error(f"创建输出文件时出错: {str(e)}")
        import traceback
        logger.error(f"详细错误: {traceback.format_exc()}")
        sys.exit(1)
    
    # 创建tabix索引
    try:
        logger.info("创建tabix索引...")
        pysam.tabix_index(output_file, preset='vcf', force=True)
        logger.info(f"索引创建完成: {output_file}.tbi")
    except Exception as e:
        logger.error(f"创建索引时出错: {str(e)}")
        # 尝试使用命令行工具
        try:
            import subprocess
            subprocess.run(['tabix', '-f', '-p', 'vcf', output_file], check=True)
            logger.info("使用命令行工具创建索引成功")
        except:
            logger.warning("无法创建tabix索引，请手动运行: tabix -f -p vcf " + output_file)
    
    # 输出统计信息
    logger.info("=== 合并统计 ===")
    logger.info(f"总变异数: {stats['total_variants']}")
    logger.info(f"跳过变异数: {stats['skipped_variants']}")
    logger.info("按caller数量分布:")
    for count, variants in sorted(stats['by_caller_count'].items()):
        logger.info(f"  {count} 个caller: {variants} 个变异")
    logger.info("按主要caller分布:")
    for caller, count in sorted(stats['by_primary_caller'].items()):
        display_name = get_caller_display_name(caller)
        logger.info(f"  {display_name}: {count} 个变异")
    
    # 验证输出文件
    try:
        with pysam.VariantFile(output_file) as check_vcf:
            output_variant_count = 0
            sample_check = None
            info_fields_found = set()
            format_fields_found = set()
            first_record_debug = True
            
            for record in check_vcf.fetch():
                output_variant_count += 1
                
                # 检查INFO字段
                for key in record.info:
                    info_fields_found.add(key)
                
                # 检查FORMAT字段和样本
                if sample_check is None and record.samples:
                    sample_names = list(record.samples.keys())
                    sample_check = sample_names
                    for sample_name in sample_names:
                        for format_key in record.samples[sample_name]:
                            format_fields_found.add(format_key)
                
                # 显示第一个记录的详细信息用于调试
                if first_record_debug:
                    logger.debug(f"第一个记录详情:")
                    logger.debug(f"  位置: {record.contig}:{record.start+1}")
                    logger.debug(f"  信息: {dict(record.info)}")
                    logger.debug(f"  样本: {list(record.samples.keys())}")
                    for sample_name in record.samples:
                        sample_data = dict(record.samples[sample_name])
                        logger.debug(f"    {sample_name}: {sample_data}")
                    first_record_debug = False
                
                # 只显示前几个记录避免过多日志
                if output_variant_count <= 3:
                    logger.debug(f"记录 {output_variant_count}: {record.contig}:{record.start+1} "
                               f"{record.ref}->{record.alts}")
            
            logger.info(f"验证输出文件: {output_variant_count} 变异")
            logger.info(f"发现INFO字段: {sorted(info_fields_found)}")
            logger.info(f"发现FORMAT字段: {sorted(format_fields_found)}")
            
            if sample_check:
                logger.info(f"样本检查通过: {sample_check}")
                
    except Exception as e:
        logger.warning(f"验证输出文件时出错: {e}")
    
    logger.info(f"合并完成，输出文件: {output_file}")

def main():
    """主函数 - 命令行模式"""
    if len(sys.argv) < 4:
        print("用法: python merge_variants.py <input1> <input2> ... <output_file> <sample_id>")
        sys.exit(1)
    
    # 解析输入参数
    input_files = sys.argv[1:-2]
    output_file = sys.argv[-2]
    sample_id = sys.argv[-1]
    
    merge_variants(input_files, output_file, sample_id)

if __name__ == "__main__":
    if 'snakemake' in globals():
        # Snakemake模式
        merge_variants(snakemake.input, snakemake.output.vcf, snakemake.params.sample_id)
    else:
        # 命令行模式
        main()