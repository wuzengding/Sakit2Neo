##scripts/Splice_Event_Specific.py
import pandas as pd
import argparse
import re
import os
import anndata as ad
from scipy.sparse import csr_matrix, csc_matrix  # 用于判断稀疏矩阵类型
from scipy import stats
from statsmodels.stats.multitest import multipletests
import numpy as np
import glob
import sys # Added for safety
import openpyxl  # 新增：用于写入 Excel

def alt_splice_TongJi_Tumor(counts_data_alt,ad_cancer):
    columns_to_keep = counts_data_alt.columns[counts_data_alt.columns != 'Alt_Splice']
    counts_data_alt = counts_data_alt.copy()
    counts_data_alt["Cancer_Mean"] = counts_data_alt[columns_to_keep].mean(axis=1)
    counts_data_alt["Cancer_Median"] = counts_data_alt[columns_to_keep].median(axis=1)
    counts_data_alt["Cancer_Min"] = counts_data_alt[columns_to_keep].min(axis=1)
    counts_data_alt["Cancer_Max"] = counts_data_alt[columns_to_keep].max(axis=1)
    # Note: If columns_to_keep is length 1 (single sample), sum is 0 or 1.
    counts_data_alt['Cancer_SampleNumber'] = (counts_data_alt[columns_to_keep] >= ad_cancer).sum(axis=1)
    counts_data_alt['Cancer_SampleRatio'] = (counts_data_alt[columns_to_keep] >= ad_cancer).sum(axis=1)/len(columns_to_keep)
    counts_data_TJ = counts_data_alt[["Alt_Splice","Cancer_Mean","Cancer_Median","Cancer_Min","Cancer_Max","Cancer_SampleNumber","Cancer_SampleRatio"]]
    return counts_data_TJ,columns_to_keep

def alt_splice_TongJi_IR_Paired(alt_splice_info,counts_data_alt,ad_cancer):
    alt_splice_info["Cancer_SampleRatio_IR_Paired"] = 0.0
    for i in alt_splice_info.index:
        alt_splice_paired = [alt_splice_info.loc[i,"Alt_Splice_x"],alt_splice_info.loc[i,"Alt_Splice_y"]]
        counts_data_alt_paired = counts_data_alt[counts_data_alt["Alt_Splice"].isin(alt_splice_paired)]
        counts_data_alt_paired = counts_data_alt_paired.drop(columns=["Alt_Splice"])
        all_positive_cols = counts_data_alt_paired.loc[:, (counts_data_alt_paired >= ad_cancer).all()]
        if len(all_positive_cols) > 0:
            pos_num = all_positive_cols.shape[1]
        else:
            pos_num = 0
        # 修复：防止除以零
        if counts_data_alt_paired.shape[1] == 0:
            IR_Paired_SampleRatio = 0.0
        else:
            IR_Paired_SampleRatio = round(pos_num/counts_data_alt_paired.shape[1],4)
        alt_splice_info.loc[i,"Cancer_SampleRatio_IR_Paired"] = IR_Paired_SampleRatio
    #print(alt_splice_info)
    alt_splice_info_x = alt_splice_info[["Alt_Splice_x","Cancer_SampleRatio_IR_Paired"]]
    alt_splice_info_y = alt_splice_info[["Alt_Splice_y","Cancer_SampleRatio_IR_Paired"]]
    alt_splice_info_x.columns = ["Alt_Splice","Cancer_SampleRatio_IR_Paired"]
    alt_splice_info_y.columns = ["Alt_Splice","Cancer_SampleRatio_IR_Paired"]
    alt_splice_info = pd.concat([alt_splice_info_x,alt_splice_info_y],axis=0)
    alt_splice_info = alt_splice_info.groupby("Alt_Splice").mean().reset_index()
    #print(alt_splice_info)
    return alt_splice_info

def read_h5ad_file(control_db):
    ctrl_db = ad.read_h5ad(os.path.join(args.database,'controls',control_db))
    ctrl_db_df = h5ad_to_dataframe(ctrl_db)
    return ctrl_db_df

def alt_splice_TongJi_Normal(ctrl_db_df):
    columns_names = ctrl_db_df.keys()
    if "Alt_Splice" in columns_names:
        columns_names = columns_names.drop("Alt_Splice")
    ctrl_db_df = ctrl_db_df.copy()
    ctrl_db_df["Mean"] = ctrl_db_df[columns_names].mean(axis=1)
    ctrl_db_df["Median"] = ctrl_db_df[columns_names].median(axis=1)
    ctrl_db_df["Min"] = ctrl_db_df[columns_names].min(axis=1)
    ctrl_db_df["Max"] = ctrl_db_df[columns_names].max(axis=1)
    ctrl_db_df['Number'] = (ctrl_db_df[columns_names] > 0).sum(axis=1)
    ctrl_db_df['Ratio'] = (ctrl_db_df[columns_names] > 0).sum(axis=1)/len(columns_names)
    return ctrl_db_df,columns_names

def alt_splice_TongJi_CSSWNormal(ctrl_db_df,ad_normal):
    columns_names = ctrl_db_df.keys()
    if "Alt_Splice" in columns_names:
        columns_names = columns_names.drop("Alt_Splice")
    ctrl_db_df = ctrl_db_df.copy()
    ctrl_db_df["Mean"] = ctrl_db_df[columns_names].mean(axis=1)
    ctrl_db_df["Median"] = ctrl_db_df[columns_names].median(axis=1)
    ctrl_db_df["Min"] = ctrl_db_df[columns_names].min(axis=1)
    ctrl_db_df["Max"] = ctrl_db_df[columns_names].max(axis=1)
    ctrl_db_df['Number'] = (ctrl_db_df[columns_names] >= ad_normal).sum(axis=1)
    ctrl_db_df['Ratio'] = (ctrl_db_df[columns_names] >= ad_normal).sum(axis=1)/len(columns_names)
    return ctrl_db_df,columns_names

def h5ad_to_dataframe(adata):
    #"""将anndata对象的X矩阵转换为DataFrame（包含行/列索引）"""
    # 1. 处理X矩阵（稀疏矩阵转密集数组）
    if isinstance(adata.X, (csr_matrix, csc_matrix)):
        X_array = adata.X.toarray()  # 稀疏矩阵转numpy数组
    else:
        X_array = adata.X  # 直接使用原数组（如已是numpy数组）
    # 2. 转换为DataFrame（行索引：样本名，列索引：特征名）
    df = pd.DataFrame(
        data=X_array,
        index=adata.obs.index,       # 行索引：样本ID（如TCGA样本名）
        columns=adata.var.index      # 列索引：特征名（如junction名称）
    )
    return df

def diff_Test(test_data,test_col_name,ctrl_data,ctrl_col_name):
    Alt_Splice = test_data["Alt_Splice"].values[0]
    test_values = list(test_data[test_col_name].values[0])
    try:
        if "Alt_Splice" in ctrl_data.keys():
            ctrl_data = ctrl_data[ctrl_data["Alt_Splice"]==Alt_Splice]
        else:
            ctrl_data = ctrl_data.loc[[Alt_Splice]]
        ctrl_mean = ctrl_data["Mean"].values[0]
        ctrl_Median = ctrl_data["Median"].values[0]
        ctrl_Min = ctrl_data["Min"].values[0]
        ctrl_Max = ctrl_data["Max"].values[0]
        ctrl_Ratio = ctrl_data["Ratio"].values[0]
        ctrl_Number = ctrl_data["Number"].values[0]
        ctrl_values = list(ctrl_data[ctrl_col_name].values[0])
    except:
        ctrl_values = []
    if len(ctrl_values) == 0:
        p_value = 0
        ctrl_mean = 0
        ctrl_Median = 0
        ctrl_Min = 0
        ctrl_Max = 0
        ctrl_Ratio = 0
        ctrl_Number = 0
    else:
        # Note: If test_values has only 1 sample, variance is 0 or undefined, t-test might warn but return result or nan.
        # For single sample, this is strictly a comparison of 1 value vs distribution.
        t, p = stats.ttest_ind(
                test_values, 
                ctrl_values,
                equal_var=False  # 默认使用Welch's t-test（方差不齐）
            )
        p_value = p
    #print(test_values)
    #print(ctrl_values)
    #print(f"{Alt_Splice}\t{len(test_values)}\t{len(ctrl_values)}\t{p_value}")
    return p_value,ctrl_mean,ctrl_Median,ctrl_Min,ctrl_Max,ctrl_Ratio,ctrl_Number 

def diff_Test_2(test_data,test_col_name,ctrl_data,ctrl_col_name):
    Alt_Splice = test_data["Alt_Splice"].values[0]
    test_values = list(test_data[test_col_name].values[0])
    try:
        ctrl_data = ctrl_data.loc[[Alt_Splice]]
        ctrl_values = list(ctrl_data[ctrl_col_name].values[0])
    except:
        ctrl_values = []
    if len(ctrl_values) == 0:
        p_value = 0
    else:
        t, p = stats.ttest_ind(
                test_values, 
                ctrl_values,
                equal_var=False  # 默认使用Welch's t-test（方差不齐）
            )
        p_value = p
    #print(f"{Alt_Splice}\t{len(test_values)}\t{len(ctrl_values)}\t{p_value}")
    #sys.exit()
    return p_value  

def calculate_signature(purned_counts_data_TJ,purned_counts_data,columns_to_keep,gtex_ctrl_df,gtex_columns_names,gtex_skin_ctrl_df,gtex_skin_columns_names,tcga_ctrl_df,tcga_columns_names,Cancer_Data,cssw_ctrl_alt_TJ,cssw_columns_names):
    purned_counts_data_TJ = purned_counts_data_TJ.copy()
    #purned_counts_data_TJ = purned_counts_data_TJ.head(1)
    all_control_number = len(gtex_columns_names) + len(gtex_skin_columns_names) + len(tcga_columns_names)
    purned_counts_data_TJ["Cancer_alt_dataset_num"] = 0
    purned_counts_data_TJ["Cancer_all_dataset_num"] = 0
    purned_counts_data_TJ["Cancer_dataset_ratio"] = 0.0
    purned_counts_data_TJ["Skin_Mean"] = 0.0
    purned_counts_data_TJ["Skin_Median"] = 0
    purned_counts_data_TJ["Skin_Min"] = 0
    purned_counts_data_TJ["Skin_Max"] = 0
    purned_counts_data_TJ["Skin_SampleNumber"] = 0
    purned_counts_data_TJ["Skin_SampleRatio"] = 0.0
    purned_counts_data_TJ["Skin_Pvalue"] = 0.0
    purned_counts_data_TJ["Skin_adjPvalue"] = 0.0
    purned_counts_data_TJ["Gtex_Mean"] = 0.0
    purned_counts_data_TJ["Gtex_Median"] = 0
    purned_counts_data_TJ["Gtex_Min"] = 0
    purned_counts_data_TJ["Gtex_Max"] = 0
    purned_counts_data_TJ["Gtex_SampleNumber"] = 0
    purned_counts_data_TJ["Gtex_SampleRatio"] = 0.0
    purned_counts_data_TJ["Gtex_Pvalue"] = 0.0
    purned_counts_data_TJ["Gtex_adjPvalue"] = 0.0
    purned_counts_data_TJ["TCGA_Mean"] = 0.0
    purned_counts_data_TJ["TCGA_Median"] = 0
    purned_counts_data_TJ["TCGA_Min"] = 0
    purned_counts_data_TJ["TCGA_Max"] = 0
    purned_counts_data_TJ["TCGA_SampleNumber"] = 0
    purned_counts_data_TJ["TCGA_SampleRatio"] = 0.0
    purned_counts_data_TJ["TCGA_Pvalue"] = 0.0
    purned_counts_data_TJ["TCGA_adjPvalue"] = 0.0
    purned_counts_data_TJ["All_Normal_SampleRatio"] = 0.0
    purned_counts_data_TJ["All_Pvalue"] = 0.0
    purned_counts_data_TJ["All_adjPvalue"] = 0.0
    purned_counts_data_TJ["Diff_Gtex"] = 0.0
    purned_counts_data_TJ["Diff_TCGA"] = 0.0
    purned_counts_data_TJ["Diff_Skin"] = 0.0
    purned_counts_data_TJ["All_CSSWNormal_SamplePosNum"] = 0.0
    purned_counts_data_TJ["All_CSSWNormal_SampleRatio"] = 0.0

    #print(purned_counts_data_TJ)
    for i in purned_counts_data_TJ.index:
        Alt_Splice = purned_counts_data_TJ.loc[i,"Alt_Splice"]
        tumor_data = purned_counts_data[purned_counts_data["Alt_Splice"]==Alt_Splice]
        #Alt_Splice = tumor_data["Alt_Splice"].values[0]
        ######
        #print(gtex_ctrl_df.loc[[Alt_Splice]])
        ######
        p_value_gtex_ctrl,ctrl_mean,ctrl_Median,ctrl_Min,ctrl_Max,ctrl_Ratio,ctrl_Number = diff_Test(tumor_data,columns_to_keep,gtex_ctrl_df,gtex_columns_names)
        purned_counts_data_TJ.loc[i,"Gtex_Pvalue"] = p_value_gtex_ctrl
        purned_counts_data_TJ.loc[i,"Gtex_Mean"] = ctrl_mean
        purned_counts_data_TJ.loc[i,"Gtex_Median"] = ctrl_Median
        purned_counts_data_TJ.loc[i,"Gtex_Min"] = ctrl_Min
        purned_counts_data_TJ.loc[i,"Gtex_Max"] = ctrl_Max
        purned_counts_data_TJ.loc[i,"Gtex_SampleRatio"] = ctrl_Ratio
        purned_counts_data_TJ.loc[i,"Gtex_SampleNumber"] = ctrl_Number
        #print(ctrl_mean)
        ######
        p_value_gtex_skin,ctrl_mean,ctrl_Median,ctrl_Min,ctrl_Max,ctrl_Ratio,ctrl_Number = diff_Test(tumor_data,columns_to_keep,gtex_skin_ctrl_df,gtex_skin_columns_names)
        purned_counts_data_TJ.loc[i,"Skin_Pvalue"] = p_value_gtex_skin
        purned_counts_data_TJ.loc[i,"Skin_Mean"] = ctrl_mean
        purned_counts_data_TJ.loc[i,"Skin_Median"] = ctrl_Median
        purned_counts_data_TJ.loc[i,"Skin_Min"] = ctrl_Min
        purned_counts_data_TJ.loc[i,"Skin_Max"] = ctrl_Max
        purned_counts_data_TJ.loc[i,"Skin_SampleRatio"] = ctrl_Ratio
        purned_counts_data_TJ.loc[i,"Skin_SampleNumber"] = ctrl_Number
        ######
        p_value_tcga,ctrl_mean,ctrl_Median,ctrl_Min,ctrl_Max,ctrl_Ratio,ctrl_Number = diff_Test(tumor_data,columns_to_keep,tcga_ctrl_df,tcga_columns_names)
        purned_counts_data_TJ.loc[i,"TCGA_Pvalue"] = p_value_tcga
        purned_counts_data_TJ.loc[i,"TCGA_Mean"] = ctrl_mean
        purned_counts_data_TJ.loc[i,"TCGA_Median"] = ctrl_Median
        purned_counts_data_TJ.loc[i,"TCGA_Min"] = ctrl_Min
        purned_counts_data_TJ.loc[i,"TCGA_Max"] = ctrl_Max
        purned_counts_data_TJ.loc[i,"TCGA_SampleRatio"] = ctrl_Ratio
        purned_counts_data_TJ.loc[i,"TCGA_SampleNumber"] = ctrl_Number
        #########################################################
        #########################################################
        p_value_cssw,ctrl_mean,ctrl_Median,ctrl_Min,ctrl_Max,ctrl_Ratio,ctrl_Number = diff_Test(tumor_data,columns_to_keep,cssw_ctrl_alt_TJ,cssw_columns_names)
        purned_counts_data_TJ.loc[i,"All_CSSWNormal_SamplePosNum"] = ctrl_Number
        purned_counts_data_TJ.loc[i,"All_CSSWNormal_SampleRatio"] = ctrl_Ratio
        #########################################################
        #########################################################
        All_Sub = []
        All_Col_Names = []
        try:
            gtex_ctrl_df_sub = gtex_ctrl_df.loc[[Alt_Splice]]
        except:
            gtex_ctrl_df_sub = []
            #gtex_ctrl_df_sub = pd.DataFrame([[0]*len(gtex_columns_names)], index=[Alt_Splice], columns=gtex_columns_names)
        try:
            gtex_skin_ctrl_df_sub = gtex_skin_ctrl_df.loc[[Alt_Splice]]
        except:
            gtex_skin_ctrl_df_sub = []
            #gtex_skin_ctrl_df_sub = pd.DataFrame([[0]*len(gtex_skin_columns_names)], index=[Alt_Splice], columns=gtex_skin_columns_names)
        try:
            tcga_ctrl_df_sub = tcga_ctrl_df.loc[[Alt_Splice]]
        except:
            tcga_ctrl_df_sub = []
            #tcga_ctrl_df_sub = pd.DataFrame([[0]*len(tcga_columns_names)], index=[Alt_Splice], columns=tcga_columns_names)
        if len(gtex_ctrl_df_sub) > 0:
            All_Col_Names = All_Col_Names + list(gtex_columns_names)
            if len(All_Sub) > 0:
                All_Sub = pd.concat([All_Sub, gtex_ctrl_df_sub], axis=1)
            else:
                All_Sub = gtex_ctrl_df_sub
        if len(gtex_skin_ctrl_df_sub) > 0:
            All_Col_Names = All_Col_Names + list(gtex_skin_columns_names)
            if len(All_Sub) > 0:
                All_Sub = pd.concat([All_Sub, gtex_skin_ctrl_df_sub], axis=1)
            else:
                All_Sub = gtex_skin_ctrl_df_sub
        if len(tcga_ctrl_df_sub) > 0:
            All_Col_Names = All_Col_Names + list(tcga_columns_names)
            if len(All_Sub) > 0:
                All_Sub = pd.concat([All_Sub, tcga_ctrl_df_sub], axis=1)
            else:
                All_Sub = tcga_ctrl_df_sub
        #print("##############################################")
        #print(All_Sub)
        if len(All_Sub) > 0:
            All_Normal_SampleRatio = ((All_Sub[All_Col_Names] > 0).sum(axis=1)/all_control_number).values[0]
        else:
            All_Normal_SampleRatio = 0.0
        #print(All_Normal_SampleRatio)
        p_value_all = diff_Test_2(tumor_data,columns_to_keep,All_Sub,All_Col_Names)
        purned_counts_data_TJ.loc[i,"All_Pvalue"] = p_value_all
        #########################################################
        #########################################################
        #########################################################
        #print(tumor_data[columns_to_keep])
        #print(tumor_data[columns_to_keep] > 0)
        mask = tumor_data[columns_to_keep] > 0
        
        # 修复：处理单行/单列情况
        if mask.empty or mask.shape[0] == 0:
            colinfo_List = columns_to_keep.tolist() if hasattr(columns_to_keep, 'tolist') else list(columns_to_keep)
        else:
            # 安全地获取满足条件的列名
            try:
                # 方法1：使用 apply 并检查结果
                colinfo = mask.apply(lambda row: row.index[row].tolist() if row.any() else [], axis=1)
                colinfo_values = colinfo.values
                if len(colinfo_values) > 0 and colinfo_values[0]:
                    colinfo_List = colinfo_values[0]
                else:
                    colinfo_List = columns_to_keep.tolist() if hasattr(columns_to_keep, 'tolist') else list(columns_to_keep)
            except (ValueError, IndexError):
                # 方法2：回退到简单逻辑
                colinfo_List = columns_to_keep.tolist() if hasattr(columns_to_keep, 'tolist') else list(columns_to_keep)
        # 修复结束
        
        # NOTE: For single sample, 'columns_to_keep' contains 1 sample. 
        # 'Cancer_Data' must contain this sample ID.
        all_dataset_info = Cancer_Data[Cancer_Data["sample_bed_id"].isin(columns_to_keep)]
        alt_dataset_info = Cancer_Data[Cancer_Data["sample_bed_id"].isin(colinfo_List)]
        
        # Safety check if Cancer_Data is missing info for the sample
        if len(all_dataset_info) == 0:
            all_dataset_num = 1
        else:
            all_dataset_num = len(np.unique(all_dataset_info["BioProject"].values))
            
        if len(alt_dataset_info) == 0:
            alt_dataset_num = 1
        else:
            alt_dataset_num = len(np.unique(alt_dataset_info["BioProject"].values))
            
        dataset_ratio = round(alt_dataset_num/all_dataset_num,4)
        purned_counts_data_TJ.loc[i,"Cancer_alt_dataset_num"] = alt_dataset_num
        purned_counts_data_TJ.loc[i,"Cancer_all_dataset_num"] = all_dataset_num
        purned_counts_data_TJ.loc[i,"Cancer_dataset_ratio"] = dataset_ratio
    ####################################################
    purned_counts_data_TJ["Diff_Gtex"] = purned_counts_data_TJ["Cancer_Mean"] - purned_counts_data_TJ["Gtex_Mean"]
    purned_counts_data_TJ["Diff_TCGA"] = purned_counts_data_TJ["Cancer_Mean"] - purned_counts_data_TJ["TCGA_Mean"]
    purned_counts_data_TJ["Diff_Skin"] = purned_counts_data_TJ["Cancer_Mean"] - purned_counts_data_TJ["Skin_Mean"]
    ####################################################
    try:
        _,Gtex_Pvalue_List_Adj,_,_ = multipletests(purned_counts_data_TJ["Gtex_Pvalue"].values,alpha=0.05,method='fdr_bh')
        _,Skin_Pvalue_List_Adj,_,_ = multipletests(purned_counts_data_TJ["Skin_Pvalue"].values,alpha=0.05,method='fdr_bh')
        _,TCGA_Pvalue_List_Adj,_,_ = multipletests(purned_counts_data_TJ["TCGA_Pvalue"].values,alpha=0.05,method='fdr_bh')
        _,All_Pvalue_List_Adj,_,_ = multipletests(purned_counts_data_TJ["All_Pvalue"].values,alpha=0.05,method='fdr_bh')
        purned_counts_data_TJ["Skin_adjPvalue"] = Skin_Pvalue_List_Adj
        purned_counts_data_TJ["Gtex_adjPvalue"] = Gtex_Pvalue_List_Adj
        purned_counts_data_TJ["TCGA_adjPvalue"] = TCGA_Pvalue_List_Adj
        purned_counts_data_TJ["All_adjPvalue"] = All_Pvalue_List_Adj
    except:
        purned_counts_data_TJ["Skin_adjPvalue"] = purned_counts_data_TJ["Skin_Pvalue"]
        purned_counts_data_TJ["Gtex_adjPvalue"] = purned_counts_data_TJ["Gtex_Pvalue"]
        purned_counts_data_TJ["TCGA_adjPvalue"] = purned_counts_data_TJ["TCGA_Pvalue"]
        purned_counts_data_TJ["All_adjPvalue"] = purned_counts_data_TJ["All_Pvalue"]
    return purned_counts_data_TJ

def IR_Specific_Filter(splice_data):
    high_specific_cutoff = 0.005
    high_specific_cutoff_cssw = 0.01
    Cancer_SampleRatio_cutoff = 0.02
    
    # MODIFICATION FOR SINGLE SAMPLE:
    # Changed recurrence from >= 2 to >= 1
    splice_data = splice_data[splice_data["Cancer_alt_dataset_num_x"]>=1]
    splice_data = splice_data[splice_data["Cancer_alt_dataset_num_y"]>=1]
    
    # For single sample, ratio is likely 1.0 (if present) or 0.0. 1.0 >= 0.02 passes.
    splice_data = splice_data[splice_data["Cancer_SampleRatio_x"]>=Cancer_SampleRatio_cutoff]
    splice_data = splice_data[splice_data["Cancer_SampleRatio_y"]>=Cancer_SampleRatio_cutoff]
    ########################################################################
    splice_data_specific = splice_data[splice_data["Diff_Gtex_x"]>0]
    splice_data_specific = splice_data_specific[splice_data_specific["Diff_TCGA_x"]>0]
    splice_data_specific = splice_data_specific[splice_data_specific["Diff_Skin_x"]>0]
    splice_data_specific = splice_data_specific[splice_data_specific["Diff_Gtex_y"]>0]
    splice_data_specific = splice_data_specific[splice_data_specific["Diff_TCGA_y"]>0]
    splice_data_specific = splice_data_specific[splice_data_specific["Diff_Skin_y"]>0]
    #splice_data_specific = splice_data_specific[splice_data_specific["All_adjPvalue_x"]<0.05]
    #splice_data_specific = splice_data_specific[splice_data_specific["All_adjPvalue_y"]<0.05]
    ###################################################################
    splice_data_high_specific = splice_data_specific[splice_data_specific["Skin_SampleRatio_x"]<high_specific_cutoff]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["Gtex_SampleRatio_x"]<high_specific_cutoff]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["TCGA_SampleRatio_x"]<high_specific_cutoff]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["Skin_SampleRatio_y"]<high_specific_cutoff]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["Gtex_SampleRatio_y"]<high_specific_cutoff]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["TCGA_SampleRatio_y"]<high_specific_cutoff]
    ###################################################################
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["All_CSSWNormal_SampleRatio_x"]<high_specific_cutoff_cssw]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["All_CSSWNormal_SampleRatio_y"]<high_specific_cutoff_cssw]
    return splice_data,splice_data_specific,splice_data_high_specific

def Alt_Specific_Filter(splice_data):
    high_specific_cutoff = 0.005
    high_specific_cutoff_cssw = 0.01
    Cancer_SampleRatio_cutoff = 0.02
    
    # MODIFICATION FOR SINGLE SAMPLE:
    # Changed recurrence from >= 2 to >= 1
    splice_data = splice_data[splice_data["Cancer_alt_dataset_num"]>=1]
    
    splice_data = splice_data[splice_data["Cancer_SampleRatio"]>=Cancer_SampleRatio_cutoff]
    ########################################################################
    splice_data_specific = splice_data[splice_data["Diff_Gtex"]>0]
    splice_data_specific = splice_data_specific[splice_data_specific["Diff_TCGA"]>0]
    splice_data_specific = splice_data_specific[splice_data_specific["Diff_Skin"]>0]
    #splice_data_specific = splice_data_specific[splice_data_specific["All_adjPvalue"]<0.05]
    ########################################################################
    splice_data_high_specific = splice_data_specific[splice_data_specific["Skin_SampleRatio"]<high_specific_cutoff]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["Gtex_SampleRatio"]<high_specific_cutoff]
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["TCGA_SampleRatio"]<high_specific_cutoff]
    ###################################################################
    splice_data_high_specific = splice_data_high_specific[splice_data_high_specific["All_CSSWNormal_SampleRatio"]<high_specific_cutoff_cssw]
    return splice_data,splice_data_specific,splice_data_high_specific

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--Splice_File", help="Splice_Event_File_List", required=True)
    parser.add_argument("-o", "--outputdir", help="output dir", required=True)
    parser.add_argument("-p", "--prunedcount", help="prunedcount", required=True)
    parser.add_argument("-d", "--database", help="database", required=True)
    parser.add_argument("-c", "--cssw_control", help="cssw_control", required=True)
    parser.add_argument("-u", "--info_file", help="info_file", required=True)
    parser.add_argument("-n", "--sample_name", help="sample name for output file naming", required=True)  # 新增参数
    parser.add_argument("-a", "--ad_normal", help="normal sample read support cutoff", required=False, default=3)
    parser.add_argument("-b", "--ad_cancer", help="cancer sample read support cutoff", required=False, default=3)
    args = parser.parse_args()
    ################################################################################################
    ################################################################################################
    ad_normal = int(args.ad_normal)
    ad_cancer = int(args.ad_cancer)
    # intergenic_cancer removed as it was not in argparse arguments
    ################################################################################################
    ################################################################################################
    alt_type_list = ["intron-retention","cassette-exon","alt-3'","alt-5'","alt-C-term","altPromoter","trans-splicing"]
    #alt_type_list = ["intron-retention","cassette-exon"]
    ################################################################################################
    purned_counts_data = pd.read_csv(args.prunedcount,sep="\t",header=0,skiprows=0)
    purned_counts_data = purned_counts_data.rename(columns={'Unnamed: 0': 'Alt_Splice'})
    #test = ["ENSG00000081923:I12.1-E13.1","ENSG00000081923:E12.1-I12.1","ENSG00000163872:E21.1-I21.1","ENSG00000163872:I21.1-E22.1","ENSG00000128487:E23.1-E24.1_20253507","ENSG00000196935:E3.1_63989964-E4.1","ENSG00000102780:E35.1-E37.1"]
    #purned_counts_data = purned_counts_data[purned_counts_data["Alt_Splice"].isin(test)]
    ################################################################################################
    
    # MODIFICATION: Read the Info file to define QCReport_data
    QCReport_data = pd.read_csv(args.info_file, sep="\t")
    
    QCReport_data["sample_bed_id"] = QCReport_data["SampleID"] + ".Aligned.sortedByCoord.out.filter.bed"
    sample_bed_id_no_contamination = ['Alt_Splice'] + list(QCReport_data["sample_bed_id"].values)
    
    # Ensure we only keep columns that exist in the pruned data
    valid_cols = [col for col in sample_bed_id_no_contamination if col in purned_counts_data.columns]
    purned_counts_data = purned_counts_data[valid_cols]
    
    ################################################################################################
    gtex_skin_ctrl_df = read_h5ad_file("gtex_skin_count.h5ad")
    gtex_skin_ctrl_df = gtex_skin_ctrl_df.loc[gtex_skin_ctrl_df.index.isin(list(purned_counts_data["Alt_Splice"].values))]
    gtex_ctrl_df = read_h5ad_file("GTEx_junction_counts.h5ad")
    gtex_ctrl_df = gtex_ctrl_df.loc[gtex_ctrl_df.index.isin(list(purned_counts_data["Alt_Splice"].values))]
    tcga_ctrl_df = read_h5ad_file("tcga_matched_control_junction_count.h5ad")
    tcga_ctrl_df = tcga_ctrl_df.loc[tcga_ctrl_df.index.isin(list(purned_counts_data["Alt_Splice"].values))]
    ################################################################################################
    cssw_control = args.cssw_control
    cssw_control_file_list = glob.glob(f"{cssw_control}/counts.original.full.*.txt")
    cssw_control_df = purned_counts_data[["Alt_Splice"]]
    for cssw_file in cssw_control_file_list:
        full_counts_normal_data = pd.read_csv(cssw_file,sep="\t",header=0,skiprows=0)
        full_counts_normal_data = full_counts_normal_data.rename(columns={'Unnamed: 0': 'Alt_Splice'})
        full_counts_normal_data = full_counts_normal_data[full_counts_normal_data["Alt_Splice"].isin(list(purned_counts_data["Alt_Splice"].values))]
        cssw_control_df = pd.merge(cssw_control_df,full_counts_normal_data,how="left",on="Alt_Splice")
    cssw_control_df = cssw_control_df.fillna(0)
    ################################################################################################
    # MODIFICATION: Use the pre-loaded QCReport_data
    Cancer_Data = QCReport_data.copy() 
    
    # Cancer_Data = pd.read_csv(args.info_file,sep="\t",header=0,skiprows=0) # Redundant, already read
    # Ensure sample_bed_id exists if re-reading or copying
    if "sample_bed_id" not in Cancer_Data.columns:
         Cancer_Data["sample_bed_id"] = Cancer_Data["SampleID"] + ".Aligned.sortedByCoord.out.filter.bed"
         
    Cancer_Data = Cancer_Data[Cancer_Data["sample_bed_id"].isin(sample_bed_id_no_contamination)]
    #print(purned_counts_data)
    ################################################################################################
    ################################################################################################
    alt_Splice_Filter_Out = pd.DataFrame()
    ################################################################################################
    
    # 新增：定义所有可能的类型并检测存在的文件
    all_alt_types = ["intron-retention","cassette-exon","alt-3","alt-5","alt-C-term","altPromoter","trans-splicing"]
    existing_alt_types = []
    splice_files_dict = {}
    
    with open(args.Splice_File, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                filename = os.path.basename(line)
                for alt_type in all_alt_types:
                    if f"_{alt_type}_" in filename:
                        existing_alt_types.append(alt_type)
                        splice_files_dict[alt_type] = line
                        break
    
    existing_alt_types = list(set(existing_alt_types))
    print(f"Processing alt_types: {existing_alt_types}")
    
    # 存储所有结果用于写入 Excel
    all_results = {}
    
    for alt_type in existing_alt_types:
        alt_type_Res_Raw = []
        splice_file = splice_files_dict.get(alt_type)
        
        if not splice_file or not os.path.exists(splice_file):
            print(f"Warning: File not found for {alt_type}: {splice_file}")
            continue
            
        # Check if file exists before reading to avoid crashing on missing types
        if os.path.exists(splice_file):
            splice_file_info_list = re.split(r"/",splice_file)
            # Assumes naming convention structure from previous rules
            # We use the dummy info to validate validation, or skip this check if running single mode blindly
            # For now, we trust the file list passed in -s
            splice_data = pd.read_csv(splice_file,sep="\t",header=0,skiprows=0)
            alt_type_Res_Raw.append(splice_data)

        if len(alt_type_Res_Raw) == 0:
            continue
            
        alt_type_Res = alt_type_Res_Raw[0]
        if "Report" in alt_type_Res.columns:
            alt_type_Res = alt_type_Res[alt_type_Res["Report"]=="YES"]
        
        if len(alt_type_Res) == 0:
             continue

        if alt_type == "intron-retention":
            alt_type_Res = alt_type_Res.sort_values(by=["SpliceInfo","Alt_Splice_x","Alt_Splice_y","PrunedCount_x","PrunedCount_y"],ascending=[False,False,False,False,False])
            alt_type_Res = alt_type_Res.drop_duplicates(subset=["SpliceInfo","Alt_Splice_x","Alt_Splice_y"],ignore_index=True,keep='first')
            alt_splice_info = alt_type_Res[["SpliceInfo","Alt_Splice_x","Alt_Splice_y"]]
            alt_splice_info = alt_splice_info.drop_duplicates(ignore_index=True)
            alt_splice_info = alt_splice_info.sort_values(by=["SpliceInfo"])
            purned_counts_data_alt = purned_counts_data[(purned_counts_data["Alt_Splice"].isin(alt_splice_info["Alt_Splice_x"])) | (purned_counts_data["Alt_Splice"].isin(alt_splice_info["Alt_Splice_y"]))]
            #print(purned_counts_data_alt)
            gtex_ctrl_alt = gtex_ctrl_df.loc[gtex_ctrl_df.index.isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            tcga_ctrl_alt = tcga_ctrl_df.loc[tcga_ctrl_df.index.isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            gtex_skin_ctrl_alt = gtex_skin_ctrl_df.loc[gtex_skin_ctrl_df.index.isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            cssw_ctrl_alt = cssw_control_df[cssw_control_df["Alt_Splice"].isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            ########
            purned_counts_data_alt_TJ,columns_to_keep = alt_splice_TongJi_Tumor(purned_counts_data_alt,ad_cancer)
            purned_counts_data_IR_Paired_TJ = alt_splice_TongJi_IR_Paired(alt_splice_info,purned_counts_data_alt,ad_cancer)
            purned_counts_data_alt_TJ = pd.merge(purned_counts_data_alt_TJ,purned_counts_data_IR_Paired_TJ,how="left",on="Alt_Splice")
            gtex_skin_ctrl_alt_TJ,gtex_skin_columns_names = alt_splice_TongJi_Normal(gtex_skin_ctrl_alt)
            gtex_ctrl_alt_TJ,gtex_columns_names = alt_splice_TongJi_Normal(gtex_ctrl_alt)
            tcga_ctrl_alt_TJ,tcga_columns_names = alt_splice_TongJi_Normal(tcga_ctrl_alt)
            cssw_ctrl_alt_TJ,cssw_columns_names = alt_splice_TongJi_CSSWNormal(cssw_ctrl_alt,ad_normal)
            purned_counts_data_alt_TJ = calculate_signature(purned_counts_data_alt_TJ,purned_counts_data,columns_to_keep,gtex_ctrl_alt_TJ,gtex_columns_names,gtex_skin_ctrl_alt_TJ,gtex_skin_columns_names,tcga_ctrl_alt_TJ,tcga_columns_names,Cancer_Data,cssw_ctrl_alt_TJ,cssw_columns_names)
            ########
            #print(alt_splice_info)
            purned_counts_data_alt_TJ_x = purned_counts_data_alt_TJ[purned_counts_data_alt_TJ["Alt_Splice"].isin(alt_splice_info["Alt_Splice_x"])]
            purned_counts_data_alt_TJ_y = purned_counts_data_alt_TJ[purned_counts_data_alt_TJ["Alt_Splice"].isin(alt_splice_info["Alt_Splice_y"])]
            purned_counts_data_alt_TJ_x = purned_counts_data_alt_TJ_x.add_suffix('_x')
            purned_counts_data_alt_TJ_y = purned_counts_data_alt_TJ_y.add_suffix('_y')
            #print(purned_counts_data_alt_TJ)
            #print(purned_counts_data_alt_TJ_x)
            #print(purned_counts_data_alt_TJ_y)
            ########
            alt_type_Res = pd.merge(alt_type_Res,purned_counts_data_alt_TJ_x,how="left",on="Alt_Splice_x")
            alt_type_Res = pd.merge(alt_type_Res,purned_counts_data_alt_TJ_y,how="left",on="Alt_Splice_y")
            alt_type_Res_filter,alt_type_Res_specific,alt_type_Res_high_specific = IR_Specific_Filter(alt_type_Res)
            #print(alt_type_Res)
        else:
            alt_type_Res = alt_type_Res.sort_values(by=["Alt_Splice","PrunedCount"],ascending=[False,False])
            alt_type_Res = alt_type_Res.drop_duplicates(subset=["Alt_Splice"],ignore_index=True,keep='first')
            alt_splice_info = alt_type_Res[["Alt_Splice"]]
            alt_splice_info = alt_splice_info.drop_duplicates(ignore_index=True)
            alt_splice_info = alt_splice_info.sort_values(by=["Alt_Splice"])
            purned_counts_data_alt = purned_counts_data[(purned_counts_data["Alt_Splice"].isin(alt_splice_info["Alt_Splice"]))]
            gtex_ctrl_alt = gtex_ctrl_df.loc[gtex_ctrl_df.index.isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            tcga_ctrl_alt = tcga_ctrl_df.loc[tcga_ctrl_df.index.isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            gtex_skin_ctrl_alt = gtex_skin_ctrl_df.loc[gtex_skin_ctrl_df.index.isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            cssw_ctrl_alt = cssw_control_df[cssw_control_df["Alt_Splice"].isin(list(purned_counts_data_alt["Alt_Splice"].values))]
            #print(alt_type_Res)
            purned_counts_data_alt_TJ,columns_to_keep = alt_splice_TongJi_Tumor(purned_counts_data_alt,ad_cancer)
            purned_counts_data_alt_TJ["Cancer_SampleRatio_IR_Paired"] = 0
            gtex_skin_ctrl_alt_TJ,gtex_skin_columns_names = alt_splice_TongJi_Normal(gtex_skin_ctrl_alt)
            gtex_ctrl_alt_TJ,gtex_columns_names = alt_splice_TongJi_Normal(gtex_ctrl_alt)
            tcga_ctrl_alt_TJ,tcga_columns_names = alt_splice_TongJi_Normal(tcga_ctrl_alt)
            cssw_ctrl_alt_TJ,cssw_columns_names = alt_splice_TongJi_CSSWNormal(cssw_ctrl_alt,ad_normal)
            purned_counts_data_alt_TJ = calculate_signature(purned_counts_data_alt_TJ,purned_counts_data,columns_to_keep,gtex_ctrl_alt_TJ,gtex_columns_names,gtex_skin_ctrl_alt_TJ,gtex_skin_columns_names,tcga_ctrl_alt_TJ,tcga_columns_names,Cancer_Data,cssw_ctrl_alt_TJ,cssw_columns_names)
            #print(purned_counts_data_alt_TJ)
            alt_type_Res = pd.merge(alt_type_Res,purned_counts_data_alt_TJ,how="left",on="Alt_Splice")
            alt_type_Res_filter,alt_type_Res_specific,alt_type_Res_high_specific = Alt_Specific_Filter(alt_type_Res)
        ##########################################################################################################
        alt_type_Filter_TJ = pd.DataFrame({
            "Alt_Splice_Type":[alt_type],
            "AS_Total_Number":[len(alt_type_Res_Raw)],
            "AS_Filter_number":[len(alt_type_Res_filter)],
            "AS_Specific_Number":[len(alt_type_Res_specific)],
            "AS_High_Specific_Number":[len(alt_type_Res_high_specific)],
        })
        alt_Splice_Filter_Out = pd.concat([alt_Splice_Filter_Out,alt_type_Filter_TJ],axis=0)
        
        # 存储结果用于写入 Excel
        all_results[alt_type] = {
            "filter": alt_type_Res_filter,
            "specific": alt_type_Res_specific, 
            "high_specific": alt_type_Res_high_specific
        }
        
        ##########################################################################################################
        # 保存中间文件（保留原有功能）
        alt_type_Res_filter.to_csv(f"{args.outputdir}/{args.sample_name}_{alt_type}_Specific_Raw.txt",sep="\t",index=False)
        alt_type_Res_high_specific.to_csv(f"{args.outputdir}/{args.sample_name}_{alt_type}_Specific_High.txt",sep="\t",index=False)
        #alt_type_Res_filter.head(10).to_csv(f"{args.outputdir}/{alt_type}_Specific_High.txt",sep="\t",index=False)
    ##########################################################################################################
    alt_Splice_Filter_Out.to_csv(f"{args.outputdir}/{args.sample_name}_Alt_Splice_Filter_TJ.txt",sep="\t",index=False)
    
    # 新增：写入 Excel 文件（整合所有类型）
    excel_path = f"{args.outputdir}/{args.sample_name}_Specific_High.xlsx"
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        # Sheet 1: 汇总统计
        alt_Splice_Filter_Out.to_excel(writer, sheet_name="Summary", index=False)
        
        # 每个 alt_type 写入单独的 sheet
        for alt_type, data_dict in all_results.items():
            # 清理 sheet 名称（Excel 限制）
            sheet_name = alt_type.replace("'", "").replace("-", "_")[:31]  # Excel sheet 名最大 31 字符
            
            # 写入 High Specific 结果
            if not data_dict["high_specific"].empty:
                data_dict["high_specific"].to_excel(writer, sheet_name=f"{sheet_name}_High", index=False)
            
            # 也写入 Filter 和 Specific 结果
            if not data_dict["filter"].empty:
                data_dict["filter"].to_excel(writer, sheet_name=f"{sheet_name}_Filter", index=False)
            
            if not data_dict["specific"].empty:
                data_dict["specific"].to_excel(writer, sheet_name=f"{sheet_name}_Specific", index=False)
    
    print(f"Excel file saved: {excel_path}")