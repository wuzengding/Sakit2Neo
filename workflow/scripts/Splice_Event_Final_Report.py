import pandas as pd
import argparse
import os
import sys

def process_merge(specific_df, peptide_df):
    """
    将特异性结果与肽段信息合并
    """
    # 确定合并键 (Key)
    # Specific results 对于 Intron Retention 通常有 Alt_Splice_x, Alt_Splice_y
    # Peptide results 通常只有 Alt_Splice
    
    merge_key_spec = "Alt_Splice"
    if "Alt_Splice" not in specific_df.columns:
        if "Alt_Splice_x" in specific_df.columns:
            merge_key_spec = "Alt_Splice_x"
        else:
            # 如果找不到键，直接返回原表（防止报错）
            print("Warning: Could not find Alt_Splice key in specific sheet.")
            return specific_df

    # 准备 Peptide DF (只保留关键列，避免混乱)
    # 你可以根据需要添加更多列
    cols_to_use = [
        "Alt_Splice", 
        "Transcript",
        "gene_biotype", 
        "Peptide_list_merge",          # 最终合并的肽段
        "Peptide_Confidence_Highest",  # 置信度 (1-4)
        "Peptide_Software_Confidence4", # 哪些软件支持 (4个都支持的)
        "Peptide_Software_Confidence3",
        "Peptide_Software_Confidence2",
        "Peptide_Software_Confidence1",
        "TransLateSoftware_Used"       # 原始列
    ]
    
    # 确保列存在
    available_cols = [c for c in cols_to_use if c in peptide_df.columns]
    peptide_subset = peptide_df[available_cols].copy()
    
    # 合并
    # specific_df 左连接 peptide_subset
    merged = pd.merge(
        specific_df, 
        peptide_subset, 
        left_on=merge_key_spec, 
        right_on="Alt_Splice", 
        how="left",
        suffixes=('', '_peptide')
    )
    
    # 如果 keys 不同 (例如 x)，drop 掉重复的 Alt_Splice 列
    if merge_key_spec != "Alt_Splice" and "Alt_Splice" in merged.columns:
        merged = merged.drop(columns=["Alt_Splice"])

    # 填充缺失值
    fill_values = {
        "Peptide_list_merge": "-",
        "Peptide_Confidence_Highest": 0,
        "Transcript": "-"
    }
    merged.fillna(value=fill_values, inplace=True)
    
    return merged

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--specific_xlsx", help="Specific_High.xlsx Input", required=True)
    parser.add_argument("-p", "--peptide_info", help="Peptide_Info.txt Input", required=True)
    parser.add_argument("-o", "--output_xlsx", help="Final Output XLSX", required=True)
    args = parser.parse_args()

    print(f"Loading Specificity Data: {args.specific_xlsx}")
    try:
        xls = pd.ExcelFile(args.specific_xlsx, engine='openpyxl')
    except Exception as e:
        print(f"Error reading Excel: {e}")
        sys.exit(1)

    print(f"Loading Peptide Data: {args.peptide_info}")
    if os.path.exists(args.peptide_info) and os.stat(args.peptide_info).st_size > 0:
        peptide_df = pd.read_csv(args.peptide_info, sep="\t")
    else:
        print("Warning: Peptide info file is empty or missing.")
        peptide_df = pd.DataFrame(columns=["Alt_Splice"])

    # 创建输出
    with pd.ExcelWriter(args.output_xlsx, engine='openpyxl') as writer:
        
        # 1. 处理 Summary Sheet (如果有)
        if "Summary" in xls.sheet_names:
            summary_df = pd.read_excel(xls, sheet_name="Summary")
            # Summary 不需要合并肽段信息，直接写入
            summary_df.to_excel(writer, sheet_name="Summary", index=False)
        
        # 2. 处理各类型 Sheet
        for sheet_name in xls.sheet_names:
            if sheet_name == "Summary": continue
            
            print(f"Processing sheet: {sheet_name}")
            sheet_df = pd.read_excel(xls, sheet_name=sheet_name)
            
            if len(sheet_df) == 0:
                sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)
                continue
                
            # 执行合并
            final_df = process_merge(sheet_df, peptide_df)
            
            # --- 排序逻辑 ---
            # 优先显示有肽段预测结果的 (Confidence > 0)
            # 其次按特异性 Ratio 排序
            if "Peptide_Confidence_Highest" in final_df.columns and "Cancer_SampleRatio" in final_df.columns:
                final_df.sort_values(
                    by=["Peptide_Confidence_Highest", "Cancer_SampleRatio"], 
                    ascending=[False, False], 
                    inplace=True
                )
            
            # --- 列重排 (可选: 把肽段列往前放) ---
            cols = list(final_df.columns)
            priority_cols = ["Alt_Splice", "Symbol", "Peptide_list_merge", "Peptide_Confidence_Highest"]
            # 如果是 IR，可能有 Alt_Splice_x
            if "Alt_Splice" not in cols and "Alt_Splice_x" in cols:
                priority_cols[0] = "Alt_Splice_x"
            
            new_order = []
            # 先放优先级列
            for c in priority_cols:
                if c in cols: 
                    new_order.append(c)
                    cols.remove(c)
            # 再放其他列
            new_order.extend(cols)
            final_df = final_df[new_order]

            # 写入 Sheet
            final_df.to_excel(writer, sheet_name=sheet_name, index=False)
            
    print(f"Final report saved to: {args.output_xlsx}")