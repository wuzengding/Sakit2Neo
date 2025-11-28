#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filters a VCF file based on selected variants from the review Excel report.
Reads multiple sheets to capture all manually selected variants.
"""

import pandas as pd
import pysam
import argparse
import logging
import sys
import os

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # --- 修改点：将参数名修改为与 Snakemake rule 一致的 --xlsx_report_Somatic ---
    parser.add_argument('--xlsx_report_Somatic', required=True, help='Input Excel report file (manual check).')
    parser.add_argument('--input_vcf', required=True, help='The original VEP-annotated VCF file.')
    parser.add_argument('--output_vcf', required=True, help='Path for the output filtered VCF file.')
    args = parser.parse_args()

    # 检查输入文件是否存在
    if not os.path.exists(args.xlsx_report_Somatic):
        logging.error(f"Input Excel file not found: {args.xlsx_report_Somatic}")
        sys.exit(1)

    # --- Step 1: Read Excel sheets and collect selected variants ---
    logging.info(f"Reading selected variants from {args.xlsx_report_Somatic}...")
    
    # 定义需要读取的 sheet 名称列表
    # 这里的名字必须与 generate_somatic_report.py 生成的 Sheet 名字完全一致
    target_sheets = ['low scale variants', 'Large scale variants']
    
    target_variants = set()
    
    try:
        # 注意：pd.ExcelFile 需要 'openpyxl' 库。如果报错 missing optional dependency，请 pip install openpyxl
        xl = pd.ExcelFile(args.xlsx_report_Somatic)
        existing_sheets = xl.sheet_names
        logging.info(f"Found sheets: {existing_sheets}")
        
        for sheet in target_sheets:
            if sheet in existing_sheets:
                logging.info(f"Processing sheet: {sheet}")
                df = pd.read_excel(args.xlsx_report_Somatic, sheet_name=sheet)
                
                # 检查是否有 Manual_Select 列
                if 'Manual_Select' not in df.columns:
                    logging.warning(f"Column 'Manual_Select' not found in sheet '{sheet}', skipping.")
                    continue
                
                # 数据清洗：转字符串、小写、去空格
                # fillna('no') 确保空值不报错
                df['Manual_Select'] = df['Manual_Select'].fillna('no').astype(str).str.lower().str.strip()
                
                # 筛选 'yes' 的行
                selected_df = df[df['Manual_Select'] == 'yes']
                logging.info(f"  - Found {len(selected_df)} selected variants in '{sheet}'")
                
                for _, row in selected_df.iterrows():
                    # 存入元组: (CHROM, POS, REF, ALT)
                    # 强制转换类型，防止 Excel 读取为数字/字符串混淆
                    chrom = str(row['CHROM']).strip()
                    pos = int(row['POS'])
                    ref = str(row['REF']).strip()
                    alt = str(row['ALT']).strip()
                    target_variants.add((chrom, pos, ref, alt))
            else:
                logging.info(f"Sheet '{sheet}' not found in report (this is normal if no variants of this type exist).")

    except ImportError:
        logging.error("Missing dependency: 'openpyxl'. Please install it via: pip install openpyxl")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Failed to read Excel file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    logging.info(f"Total unique variants marked as 'yes': {len(target_variants)}")

    # 如果没有选中任何变异，生成一个空的 VCF (只带 Header)
    if not target_variants:
        logging.warning("No variants marked 'yes' found. Creating empty VCF.")
        with pysam.VariantFile(args.input_vcf) as vcf_in:
            with pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header) as vcf_out:
                pass 
        try:
            pysam.tabix_index(args.output_vcf, preset="vcf")
        except Exception as e:
            logging.warning(f"Could not index empty VCF (expected if empty): {e}")
        sys.exit(0)

    # --- Step 2: Iterate VCF and filter ---
    logging.info(f"Filtering {args.input_vcf}...")
    records_written = 0
    
    try:
        with pysam.VariantFile(args.input_vcf) as vcf_in:
            with pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header) as vcf_out:
                for record in vcf_in:
                    # 检查该记录的任意一个ALT是否在目标列表中
                    write_record = False
                    
                    # 简单的 CHROM 匹配优化：如果当前的 chrom 甚至不在我们的列表中，可以快速跳过吗？
                    # 由于 VCF 通常是有序的，这里简单遍历即可，无需过度优化，除非文件巨大。
                    
                    for i, alt in enumerate(record.alts):
                        # 构造 key 进行查找
                        # 注意：record.chrom 可能是 'chr1' 或 '1'，需确保与 Excel 一致
                        # 这里假设两者格式已经一致（通常都来自同一个上游流程）
                        key = (str(record.chrom), int(record.pos), str(record.ref), str(alt))
                        
                        if key in target_variants:
                            write_record = True
                            break
                    
                    if write_record:
                        vcf_out.write(record)
                        records_written += 1
    except Exception as e:
        logging.error(f"Error processing VCF: {e}")
        sys.exit(1)

    logging.info(f"Wrote {records_written} records to {args.output_vcf}.")

     #--- Step 3: Index ---
    logging.info("Indexing output VCF...")
    try:
        pysam.tabix_index(args.output_vcf, preset="vcf")
    except Exception as e:
        logging.error(f"Failed to index VCF: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()