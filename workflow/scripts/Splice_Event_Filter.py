# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import re
import os
#python /mnt/user/dingyu/01.Shared_Neoantigen/01.SNAF/01.SNAF_Snake_Pipeline/script/Splice_Event_Filter.py -i /mnt/user/dingyu/01.Shared_Neoantigen/01.SNAF/03.Test_Analysis/Osteosarcoma_RNASeq_Cancer/07.AltAnalyze/bed/SRR10423027.Aligned.sortedByCoord.out__intronJunction.bed -j /mnt/user/dingyu/01.Shared_Neoantigen/01.SNAF/03.Test_Analysis/Osteosarcoma_RNASeq_Cancer/07.AltAnalyze/bed/SRR10423027.Aligned.sortedByCoord.out__junction.bed -e /mnt/user/dingyu/01.Shared_Neoantigen/01.SNAF/03.Test_Analysis/Osteosarcoma_RNASeq_Cancer/07.AltAnalyze/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt -s SRR10423027 -o /mnt/user/dingyu/01.Shared_Neoantigen/01.SNAF/03.Test_Analysis/Osteosarcoma_RNASeq_Cancer/SRR10423028/11.Splice_Filter -g /mnt/user/dingyu/01.Shared_Neoantigen/01.SNAF/02.database/Gtf_Bed_File_Dir_Inhouse
def process_bed_coordinates(chr, exon1_start, exon2_stop, strand, exon1_len, exon2_len):
    """
    处理BED文件坐标，完全按照AltAnalyze逻辑
    参数:
        chr: 染色体名称
        exon1_start: 第一个exon起始位置
        exon2_stop: 第二个exon结束位置  
        strand: 链方向 (+/-)
        exon1_len: 第一个exon长度
        exon2_len: 第二个exon长度
    
    返回:
        tuple: (chr, exon1_stop, exon2_start) - AltAnalyze使用的坐标键
    """
    # 染色体名称标准化 - 完全按照AltAnalyze逻辑
    if 'chr' not in chr:
        chr = 'chr' + chr
    if chr == 'chrM': 
        chr = 'chrMT'  # MT是Ensembl约定
    if chr == 'M': 
        chr = 'MT'
    
    # 坐标计算 - 完全按照AltAnalyze逻辑
    if strand == '-':  # 负链处理
        if (exon1_len + exon2_len) == 0:  # Kallisto-Splice特殊情况
            exon1_stop = exon1_start
            exon2_start = exon2_stop
        else:  # 标准情况
            exon1_stop = exon1_start + exon1_len
            exon2_start = exon2_stop - exon2_len + 1
        
        # 负链关键步骤：交换exon顺序
        a = (exon1_start, exon1_stop)
        b = (exon2_start, exon2_stop)
        exon1_stop, exon1_start = b  # 将b的值赋给exon1
        exon2_stop, exon2_start = a  # 将a的值赋给exon2
        
    else:  # 正链处理
        if (exon1_len + exon2_len) == 0:  # Kallisto-Splice特殊情况
            exon1_stop = exon1_start
            exon2_start = exon2_stop
        else:  # 标准情况
            exon1_stop = exon1_start + exon1_len
            exon2_start = exon2_stop - exon2_len + 1
    return (chr, exon1_stop, exon2_start)

def bed_file_Process(sample_data):
    for i in sample_data.index:
        chr = sample_data.loc[i,0]
        exon1_start = sample_data.loc[i,1]
        exon2_stop = sample_data.loc[i,2]
        strand = sample_data.loc[i,5]
        #IntronInfo = sample_data.loc[i,"IntronInfo"]
        #SpliceSite = sample_data.loc[i,"SpliceSite"]
        lengths = sample_data.loc[i,10]
        lengths_info = re.split(',', lengths)
        exon1_len = lengths_info[0]
        exon2_len = lengths_info[1]
        exon1_len = int(exon1_len)
        exon2_len = int(exon2_len)
        exon1_start = int(exon1_start)
        exon2_stop = int(exon2_stop)
        chr, exon1_stop, exon2_start = process_bed_coordinates(chr, exon1_start, exon2_stop, strand, exon1_len, exon2_len)
        pos= f"{chr}:{exon1_stop}-{exon2_start}"
        sample_data.loc[i,"pos"] = pos
    ##########################################################################################
    ##########################################################################################
    sample_data = sample_data[[0,3,4,5,"pos"]]
    #print(sample_data)
    sample_data.columns = ["chr","junction_id","support_reads","strand","pos_AltAnalyze"]
    sample_data = sample_data.copy()
    sample_data["MergeKey"] = sample_data["junction_id"].str.split(':',expand=True)[0] + "_" + sample_data["pos_AltAnalyze"]
    sample_data[["IntronInfo","SpliceSite"]] = sample_data["junction_id"].str.split('-',expand=True)
    sample_data["SpliceSite"] = sample_data["SpliceSite"].astype(int)
    #print(sample_data)
    #print(sample_data)
    odd_sample = sample_data[sample_data.index % 2 == 0]
    even_sample = sample_data[sample_data.index % 2 == 1]
    odd_sample.columns = ["chr_x","junction_id_x","support_reads_x","strand_x","pos_AltAnalyze_x","MergeKey_x","IntronInfo_x","SpliceSite_x"]
    even_sample.columns = ["chr_y","junction_id_y","support_reads_y","strand_y","pos_AltAnalyze_y","MergeKey_y","IntronInfo_y","SpliceSite_y"]
    combined_data = pd.concat([odd_sample.reset_index(drop=True), even_sample.reset_index(drop=True)], axis=1)
    combined_data["pos_Splice"] = combined_data["chr_x"] + ":" + (combined_data[["SpliceSite_x","SpliceSite_y"]].min(axis=1)).astype(str) + "-" + (combined_data[["SpliceSite_x","SpliceSite_y"]].max(axis=1)).astype(str)
    combined_data["SpliceInfo"] = combined_data["IntronInfo_x"] + "_" + combined_data["pos_Splice"]
    #print(odd_sample)
    #print(even_sample)
    #print(combined_data)
    #combined_data_1 = combined_data[["MergeKey_x","pos_AltAnalyze_x","SpliceInfo","pos_Raw","strand_x"]]
    #combined_data_2 = combined_data[["MergeKey_y","pos_AltAnalyze_y","SpliceInfo","pos_Raw","strand_y"]]
    #combined_data_1.columns = ["MergeKey","pos_Examined","SpliceInfo","pos_Raw","Strand"]
    #combined_data_2.columns = ["MergeKey","pos_Examined","SpliceInfo","pos_Raw","Strand"]
    combined_data_1 = combined_data[["MergeKey_x","SpliceInfo","pos_Splice","strand_x","support_reads_x"]]
    combined_data_2 = combined_data[["MergeKey_y","SpliceInfo","pos_Splice","strand_y","support_reads_y"]]
    combined_data_1.columns = ["MergeKey","SpliceInfo","pos_Splice","Strand","support_reads_intronbed"]
    combined_data_2.columns = ["MergeKey","SpliceInfo","pos_Splice","Strand","support_reads_intronbed"]
    combined_data = pd.concat([combined_data_1, combined_data_2], axis=0, ignore_index= True)
    #print(combined_data_1)
    #print(combined_data_2)
    return combined_data

def bed_file_Process_2(sample_data):
    for i in sample_data.index:
        chr = sample_data.loc[i,0]
        exon1_start = sample_data.loc[i,1]
        exon2_stop = sample_data.loc[i,2]
        strand = sample_data.loc[i,5]
        #IntronInfo = sample_data.loc[i,"IntronInfo"]
        #SpliceSite = sample_data.loc[i,"SpliceSite"]
        lengths = sample_data.loc[i,10]
        lengths_info = re.split(',', lengths)
        exon1_len = lengths_info[0]
        exon2_len = lengths_info[1]
        exon1_len = int(exon1_len)
        exon2_len = int(exon2_len)
        exon1_start = int(exon1_start)
        exon2_stop = int(exon2_stop)
        chr, exon1_stop, exon2_start = process_bed_coordinates(chr, exon1_start, exon2_stop, strand, exon1_len, exon2_len)
        pos= f"{chr}:{exon1_stop}-{exon2_start}"
        sample_data.loc[i,"pos"] = pos
    ##########################################################################################
    ##########################################################################################
    sample_data = sample_data[[0,1,2,3,4,5,"pos"]]
    #print(sample_data)
    sample_data.columns = ["chr","start_min","stop_max","junction_id","support_reads_junctionbed","Strand","pos_AltAnalyze"]
    sample_data = sample_data[["pos_AltAnalyze","junction_id","Strand","support_reads_junctionbed","start_min","stop_max"]]
    sample_data.columns = ["pos_Examined","junction_id","Strand","start_min","stop_max","support_reads_junctionbed"]
    return sample_data

def IR_bed(AnnoMergeData,outdir,SampleID,bedfiledir):
    outfile1 = f"{outdir}/{SampleID}_IR_Anno_Merge_inner.bed"
    outfile2 = f"{outdir}/{SampleID}_IR_Anno_Merge_pos1.bed"
    outfile3 = f"{outdir}/{SampleID}_IR_Anno_Merge_pos2.bed"
    AS_Pair = AnnoMergeData[["SpliceInfo","Alt_Splice_x","Alt_Splice_y","Strand_x","Strand_y","pos_Examined_x","pos_Examined_y","pos_IGV"]]
    AS_Pair = AS_Pair.drop_duplicates()
    if len(AS_Pair) > 0:
        AS_Pair[['chr', 'start', 'end']] = AS_Pair['pos_IGV'].str.split(':|-', expand=True)
        AS_Pair["start_left"] = (AS_Pair["start"].astype(int)) - 1  
        AS_Pair["end_left"] = (AS_Pair["start"].astype(int)) + 1  
        AS_Pair["start_right"] = (AS_Pair["end"].astype(int)) - 1  
        AS_Pair["end_right"] = (AS_Pair["end"].astype(int)) + 1  
        AS_Pair[['chr', 'start', 'end', 'pos_IGV', "SpliceInfo"]].to_csv(outfile1,sep="\t",header=False,index=False)
        AS_Pair[['chr', 'start_left', 'end_left', 'pos_IGV', "SpliceInfo"]].to_csv(outfile2,sep="\t",header=False,index=False)
        AS_Pair[['chr', 'start_right', 'end_right', 'pos_IGV', "SpliceInfo"]].to_csv(outfile3,sep="\t",header=False,index=False)
        intersect_Script1 = f"bedtools intersect -wao -a {outfile2} -b {bedfiledir}/exon.rightpos.bed > {outdir}/{SampleID}_IR_pos1.intersect.bed"
        intersect_Script2 = f"bedtools intersect -wao -a {outfile3} -b {bedfiledir}/exon.leftpos.bed > {outdir}/{SampleID}_IR_pos2.intersect.bed"
        intersect_Script3 = f"bedtools intersect -wao -a {outfile1} -b {bedfiledir}/intron.bed > {outdir}/{SampleID}_IR_intron.intersect.bed"
        intersect_Script4 = f"bedtools intersect -wao -a {outfile1} -b {bedfiledir}/exon.bed > {outdir}/{SampleID}_IR_exon.intersect.bed"
        intersect_Script5 = f"bedtools intersect -wao -a {outfile1} -b {bedfiledir}/exon.sort.merge.bed > {outdir}/{SampleID}_IR_exonMerge.intersect.bed"
        os.system(intersect_Script1)
        os.system(intersect_Script2)
        os.system(intersect_Script3)
        os.system(intersect_Script4)
        os.system(intersect_Script5)
        IR_Filter_Res = IR_bed_Filter(f"{outdir}/{SampleID}_IR_intron.intersect.bed",
                                        f"{outdir}/{SampleID}_IR_exon.intersect.bed",
                                        f"{outdir}/{SampleID}_IR_exonMerge.intersect.bed",
                                        f"{outdir}/{SampleID}_IR_pos1.intersect.bed",
                                        f"{outdir}/{SampleID}_IR_pos2.intersect.bed")
    else:
        IR_Filter_Res = pd.DataFrame()
    return IR_Filter_Res

def IR_bed_Filter(intron_intersect_file,exon_intersect_file,exon_intersect_merge_file,pos1_intersect_file,pos2_intersect_file):
    #intron_intersect_file = f"{outdir}/05.SampleSplice_IR_intron.intersect.bed"
    intron_intersect_data = pd.read_csv(intron_intersect_file,sep="\t",header=None)
    outdataframe = intron_intersect_data[[4]]
    outdataframe = outdataframe.drop_duplicates()
    #print(intron_intersect_data)
    intron_intersect_data["len"] = intron_intersect_data[2] - intron_intersect_data[1]
    intron_intersect_data["intersect_intron"] = intron_intersect_data[10]/intron_intersect_data['len']
    intron_intersect_data["gene_Examined"] = intron_intersect_data[4].str.split(":",expand=True)[0]
    intron_intersect_data[["gene_anno","trans_anno"]] = intron_intersect_data[8].str.split(":",expand=True)
    #print(intron_intersect_data)
    intron_intersect_data_samegene = intron_intersect_data[intron_intersect_data["gene_Examined"]==intron_intersect_data["gene_anno"]]
    intron_intersect_data_samegene = intron_intersect_data_samegene[intron_intersect_data_samegene["intersect_intron"]>0.95]
    #print(intron_intersect_data_samegene)
    #intron_intersect_data_samegene = intron_intersect_data_samegene.groupby([4])["gene_anno"].agg(lambda x: '|'.join(x)).reset_index()
    intron_intersect_data_samegene = intron_intersect_data_samegene.groupby([4]).agg({
        'gene_anno': lambda x: '|'.join(x.astype(str)),
        'trans_anno': lambda x: '|'.join(x.astype(str))
    }).reset_index()
    #print(intron_intersect_data_samegene)

    #exon_intersect_file = f"{outdir}/05.SampleSplice_IR_exon.intersect.bed"
    exon_intersect_data = pd.read_csv(exon_intersect_file,sep="\t",header=None)
    print(exon_intersect_data)
    #print(exon_intersect_data)
    exon_intersect_data["len"] = exon_intersect_data[2] - exon_intersect_data[1]
    exon_intersect_data["intersect_exon"] = exon_intersect_data[10]/exon_intersect_data['len']
    exon_intersect_data["gene_Examined"] = exon_intersect_data[4].str.split(":",expand=True)[0]
    try:
        exon_intersect_data[["gene_anno","trans_anno"]] = exon_intersect_data[8].str.split(":",expand=True)
    except:
        exon_intersect_data[["gene_anno","trans_anno"]] = "."
    exon_intersect_data_samegene = exon_intersect_data[exon_intersect_data["gene_Examined"]==exon_intersect_data["gene_anno"]]
    exon_intersect_data_samegene = exon_intersect_data_samegene[exon_intersect_data_samegene["intersect_exon"]>0.05]
    #exon_intersect_data_samegene = exon_intersect_data_samegene.groupby([4])["gene_anno"].agg(lambda x: '|'.join(x)).reset_index()
    exon_intersect_data_samegene = exon_intersect_data_samegene.groupby([4]).agg({
        'gene_anno': lambda x: '|'.join(x.astype(str)),
        'trans_anno': lambda x: '|'.join(x.astype(str))
    }).reset_index()
    #print(exon_intersect_data_samegene)

    #exon_intersect_merge_file = f"{outdir}/05.SampleSplice_IR_exonMerge.intersect.bed"
    exon_intersect_merge_data = pd.read_csv(exon_intersect_merge_file,sep="\t",header=None)
    exon_intersect_merge_data["len"] = exon_intersect_merge_data[2] - exon_intersect_merge_data[1]
    exon_intersect_merge_data["intersect_exonmerge"] = exon_intersect_merge_data[8]/exon_intersect_merge_data['len']
    exon_intersect_merge_data = exon_intersect_merge_data.groupby([4])["intersect_exonmerge"].sum().reset_index()
    #print(exon_intersect_merge_data)

    #pos1_intersect_file = f"{outdir}/05.SampleSplice_IR_pos1.intersect.bed"
    pos1_intersect_data = pd.read_csv(pos1_intersect_file,sep="\t",header=None)
    pos1_intersect_data["gene_Examined"] = pos1_intersect_data[4].str.split(":",expand=True)[0]
    try:
        pos1_intersect_data[["gene_anno","trans_anno"]] = pos1_intersect_data[8].str.split(":",expand=True)
    except:
        pos1_intersect_data[["gene_anno","trans_anno"]] = "."
    pos1_intersect_data_samegene = pos1_intersect_data[pos1_intersect_data["gene_Examined"]==pos1_intersect_data["gene_anno"]]
    #pos1_intersect_data_samegene = pos1_intersect_data_samegene.groupby([4])["gene_anno"].agg(lambda x: '|'.join(x)).reset_index()
    pos1_intersect_data_samegene = pos1_intersect_data_samegene.groupby([4]).agg({
        'gene_anno': lambda x: '|'.join(x.astype(str)),
        'trans_anno': lambda x: '|'.join(x.astype(str))
    }).reset_index()
    #print(pos1_intersect_data_samegene)

    #pos2_intersect_file = f"{outdir}/05.SampleSplice_IR_pos2.intersect.bed"
    pos2_intersect_data = pd.read_csv(pos2_intersect_file,sep="\t",header=None)
    pos2_intersect_data["gene_Examined"] = pos2_intersect_data[4].str.split(":",expand=True)[0]
    try:
        pos2_intersect_data[["gene_anno","trans_anno"]] = pos2_intersect_data[8].str.split(":",expand=True)
    except:
        pos2_intersect_data[["gene_anno","trans_anno"]] = "."
    pos2_intersect_data_samegene = pos2_intersect_data[pos2_intersect_data["gene_Examined"]==pos2_intersect_data["gene_anno"]]
    #pos2_intersect_data_samegene = pos2_intersect_data_samegene.groupby([4])["gene_anno"].agg(lambda x: '|'.join(x)).reset_index()
    pos2_intersect_data_samegene = pos2_intersect_data_samegene.groupby([4]).agg({
        'gene_anno': lambda x: '|'.join(x.astype(str)),
        'trans_anno': lambda x: '|'.join(x.astype(str))
    }).reset_index()
    #print(pos2_intersect_data_samegene)

    intron_intersect_data_samegene.columns = ["SpliceInfo","IntronGene","IntronTrans"]
    exon_intersect_data_samegene.columns = ["SpliceInfo","exonGene","exonTrans"]
    pos1_intersect_data_samegene.columns = ["SpliceInfo","Pos1Gene","Pos1Trans"]
    pos2_intersect_data_samegene.columns = ["SpliceInfo","Pos2Gene","Pos2Trans"]
    exon_intersect_merge_data.columns = ["SpliceInfo","intersect_exonmerge"]
    outdataframe.columns = ["SpliceInfo"]

    outdataframe = pd.merge(outdataframe,intron_intersect_data_samegene,on="SpliceInfo",how="left")
    outdataframe = pd.merge(outdataframe,exon_intersect_data_samegene,on="SpliceInfo",how="left")
    outdataframe = pd.merge(outdataframe,pos1_intersect_data_samegene,on="SpliceInfo",how="left")
    outdataframe = pd.merge(outdataframe,pos2_intersect_data_samegene,on="SpliceInfo",how="left")
    outdataframe = pd.merge(outdataframe,exon_intersect_merge_data,on="SpliceInfo",how="left")

    outdataframe = outdataframe.fillna(0)

    #print(outdataframe)
    outdataframe["Report"] = "YES"
    outdataframe["GeneDetail"] = "YES"
    outdataframe["TranscriptDetail"] = "YES"

    for i in outdataframe.index:
        Pos1Gene = re.split(r"\|",str(outdataframe.loc[i,"Pos1Gene"]))
        Pos2Gene = re.split(r"\|",str(outdataframe.loc[i,"Pos2Gene"]))
        IntronGene = re.split(r"\|",str(outdataframe.loc[i,"IntronGene"]))
        exonGene = re.split(r"\|",str(outdataframe.loc[i,"exonGene"]))
        intersect_exonmerge = outdataframe.loc[i,"intersect_exonmerge"]
        common_elements = set(Pos1Gene) & set(Pos2Gene) & set(IntronGene)
        diff_elements = common_elements - set(exonGene)
        diff_elements = list(diff_elements)
        if "0" in diff_elements:
            diff_elements.remove("0") 
        if len(diff_elements) > 0:
            if intersect_exonmerge <= 0.05:
                Trans = "YES"
            else:
                Trans = "NO"
            GeneDetail = "|".join(diff_elements)
        else:
            Trans = "NO"
            GeneDetail = "-"
        ################################################################
        Pos1Trans = re.split(r"\|",str(outdataframe.loc[i,"Pos1Trans"]))
        Pos2Trans = re.split(r"\|",str(outdataframe.loc[i,"Pos2Trans"]))
        IntronTrans = re.split(r"\|",str(outdataframe.loc[i,"IntronTrans"]))
        exonTrans = re.split(r"\|",str(outdataframe.loc[i,"exonTrans"]))
        common_elements_trans = set(Pos1Trans) & set(Pos2Trans) & set(IntronTrans)
        diff_elements_Trans = common_elements_trans - set(exonTrans)
        diff_elements_Trans = list(diff_elements_Trans)
        if "0" in diff_elements_Trans:
            diff_elements_Trans.remove("0") 
        if Trans == "YES":
            TranscriptDetail = "|".join(diff_elements_Trans)
        else:
            TranscriptDetail = "-"
        ################################################################
        outdataframe.loc[i,"Report"] = Trans
        outdataframe.loc[i,"GeneDetail"] = GeneDetail
        outdataframe.loc[i,"TranscriptDetail"] = TranscriptDetail
    return outdataframe

def Exon_Skip_bed(cassette_data,outdir,SampleID,bedfiledir):
    cassette_out1 = f"{outdir}/{SampleID}_Exon_Skipping_leftpos.bed"
    cassette_out2 = f"{outdir}/{SampleID}_Exon_Skipping_innerpos.bed"
    cassette_out3 = f"{outdir}/{SampleID}_Exon_Skipping_rightpos.bed"
    cassette_data[["chr","start_left", "start", "Alt_Splice","pos_IGV"]].to_csv(cassette_out1,sep="\t",header=False,index=False)
    cassette_data[["chr","start", "start_right", "Alt_Splice","pos_IGV"]].to_csv(cassette_out2,sep="\t",header=False,index=False)
    cassette_data[["chr","start_right", "end", "Alt_Splice","pos_IGV"]].to_csv(cassette_out3,sep="\t",header=False,index=False)
    #############################################################
    cassette_Script1 = f"bedtools intersect -wao -a {cassette_out1} -b {bedfiledir}/exon.rightpos.bed > {outdir}/{SampleID}_Exon_Skipping_leftpos_intersect.bed"
    cassette_Script2 = f"bedtools intersect -wao -a {cassette_out3} -b {bedfiledir}/exon.leftpos.bed > {outdir}/{SampleID}_Exon_Skipping_rightpos_intersect.bed"
    cassette_Script3 = f"bedtools intersect -wao -a {cassette_out2} -b {bedfiledir}/exon.bed > {outdir}/{SampleID}_Exon_Skipping_inner_exon_intersect.bed"
    cassette_Script4 = f"bedtools intersect -wao -a {cassette_out2} -b {bedfiledir}/intron.bed > {outdir}/{SampleID}_Exon_Skipping_inner_intron_intersect.bed"
    os.system(cassette_Script1)
    os.system(cassette_Script2)
    os.system(cassette_Script3)
    os.system(cassette_Script4)
    Exon_Skip_Filter_Res = Exon_Skip_bed_Filter(f"{outdir}/{SampleID}_Exon_Skipping_leftpos_intersect.bed",
                                                f"{outdir}/{SampleID}_Exon_Skipping_rightpos_intersect.bed",
                                                f"{outdir}/{SampleID}_Exon_Skipping_inner_exon_intersect.bed",
                                                f"{outdir}/{SampleID}_Exon_Skipping_inner_intron_intersect.bed")
    return Exon_Skip_Filter_Res

def Exon_Skip_bed_Filter(leftpos_intersect_file,rightpos_intersect_file,inner_exon_intersect_file,inner_intron_intersect_file):
    leftpos_intersect_data = pd.read_csv(leftpos_intersect_file, sep="\t", header=None, skiprows=0)
    rightpos_intersect_data = pd.read_csv(rightpos_intersect_file, sep="\t", header=None, skiprows=0)
    inner_exon_intersect_data = pd.read_csv(inner_exon_intersect_file, sep="\t", header=None, skiprows=0)
    inner_intron_intersect_data = pd.read_csv(inner_intron_intersect_file, sep="\t", header=None, skiprows=0)
    ########################################
    ########################################
    outdataframe = leftpos_intersect_data[[3]]
    outdataframe = outdataframe.drop_duplicates()
    #print(leftpos_intersect_data)
    leftpos_intersect_data["gene_Examined"] = leftpos_intersect_data[3].str.split(":",expand=True)[0]
    leftpos_intersect_data[["gene_anno","trans_anno"]] = leftpos_intersect_data[8].str.split(":",expand=True)
    leftpos_intersect_data_samegene = leftpos_intersect_data[leftpos_intersect_data["gene_Examined"]==leftpos_intersect_data["gene_anno"]]
    leftpos_intersect_data_samegene = leftpos_intersect_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #print(leftpos_intersect_data_samegene)
    ########################################
    ########################################
    #print(leftpos_intersect_data)
    rightpos_intersect_data["gene_Examined"] = rightpos_intersect_data[3].str.split(":",expand=True)[0]
    rightpos_intersect_data[["gene_anno","trans_anno"]] = rightpos_intersect_data[8].str.split(":",expand=True)
    rightpos_intersect_data_samegene = rightpos_intersect_data[rightpos_intersect_data["gene_Examined"]==rightpos_intersect_data["gene_anno"]]
    rightpos_intersect_data_samegene = rightpos_intersect_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #print(rightpos_intersect_data_samegene)
    ########################################
    ########################################
    #print(inner_exon_intersect_data)
    inner_exon_intersect_data["len"] = inner_exon_intersect_data[7] - inner_exon_intersect_data[6]
    inner_exon_intersect_data["intersect_ratio"] = inner_exon_intersect_data[10] / inner_exon_intersect_data["len"]
    inner_exon_intersect_data["gene_Examined"] = inner_exon_intersect_data[3].str.split(":",expand=True)[0]
    inner_exon_intersect_data[["gene_anno","trans_anno"]] = inner_exon_intersect_data[8].str.split(":",expand=True)
    inner_exon_intersect_data_samegene = inner_exon_intersect_data[inner_exon_intersect_data["gene_Examined"]==inner_exon_intersect_data["gene_anno"]]
    inner_exon_intersect_data_samegene = inner_exon_intersect_data_samegene[inner_exon_intersect_data_samegene["intersect_ratio"]==1]
    inner_exon_intersect_data_samegene = inner_exon_intersect_data_samegene.groupby([3,"intersect_ratio"])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #print(inner_exon_intersect_data_samegene)
    ########################################
    ########################################
    #print(inner_intron_intersect_data)
    inner_intron_intersect_data["len"] = inner_intron_intersect_data[7] - inner_intron_intersect_data[6]
    inner_intron_intersect_data["intersect_ratio"] = inner_intron_intersect_data[10] / inner_intron_intersect_data["len"]
    inner_intron_intersect_data["gene_Examined"] = inner_intron_intersect_data[3].str.split(":",expand=True)[0]
    inner_intron_intersect_data[["gene_anno","trans_anno"]] = inner_intron_intersect_data[8].str.split(":",expand=True)
    inner_intron_intersect_data_samegene = inner_intron_intersect_data[inner_intron_intersect_data["gene_Examined"]==inner_intron_intersect_data["gene_anno"]]
    inner_intron_intersect_data_samegene = inner_intron_intersect_data_samegene[inner_intron_intersect_data_samegene["intersect_ratio"]==1]
    inner_intron_intersect_data_samegene = inner_intron_intersect_data_samegene.groupby([3,"trans_anno"]).agg({"intersect_ratio":"sum"}).reset_index()
    inner_intron_intersect_data_samegene = inner_intron_intersect_data_samegene[inner_intron_intersect_data_samegene["intersect_ratio"]>=2]
    inner_intron_intersect_data_samegene = inner_intron_intersect_data_samegene.groupby([3])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    inner_intron_intersect_data_samegene["intersect_ratio"] = 2
    #print(inner_intron_intersect_data_samegene)
    #print(inner_intron_intersect_data_diffgene)
    ########################################
    ########################################
    outdataframe = pd.merge(outdataframe,leftpos_intersect_data_samegene,how='left',on=[3])
    outdataframe = pd.merge(outdataframe,rightpos_intersect_data_samegene,how='left',on=[3])
    outdataframe.columns = [3,"leftpos_samegene","leftpos_samegene_trans","rightpos_samegene","rightpos_samegene_trans"]
    outdataframe = pd.merge(outdataframe,inner_exon_intersect_data_samegene,how='left',on=[3])
    outdataframe = pd.merge(outdataframe,inner_intron_intersect_data_samegene,how='left',on=[3])
    outdataframe.columns = ["Alt_Splice","leftpos_samegene","leftpos_samegene_trans","rightpos_samegene","rightpos_samegene_trans","inner_exon_samegene","inner_exon_samegene_trans","inner_intron_samegene_trans","inner_intron_samegene"]
    outdataframe = outdataframe.fillna(0)
    outdataframe["Report"] = "YES"
    outdataframe["TranscriptDetail"] = "YES"
    for i in outdataframe.index:
        leftpos_samegene_trans = re.split(r"\|",str(outdataframe.loc[i,"leftpos_samegene_trans"]))
        rightpos_samegene_trans = re.split(r"\|",str(outdataframe.loc[i,"rightpos_samegene_trans"]))
        inner_exon_samegene_trans = re.split(r"\|",str(outdataframe.loc[i,"inner_exon_samegene_trans"]))
        inner_intron_samegene_trans = re.split(r"\|",str(outdataframe.loc[i,"inner_intron_samegene_trans"]))
        common_elements = set(leftpos_samegene_trans) & set(rightpos_samegene_trans) & set(inner_exon_samegene_trans) & set(inner_intron_samegene_trans)
        common_elements = list(common_elements)
        if "0" in common_elements:
            common_elements.remove("0") 
        if len(common_elements) > 0:
            Trans = "YES"
            TransDetail = "|".join(common_elements)
        else:
            Trans = "NO"
            TransDetail = "-"
        outdataframe.loc[i,"Report"] = Trans
        outdataframe.loc[i,"TranscriptDetail"] = TransDetail
    #print(outdataframe)
    return outdataframe

def Alt_3_bed(ID_Data_3,outdir,SampleID,bedfiledir):
    ID_Data_3_out1 = f"{outdir}/{SampleID}_Alt-3-5pos.bed"
    ID_Data_3_out2 = f"{outdir}/{SampleID}_Alt-3-3pos.bed"
    ID_Data_3_out3 = f"{outdir}/{SampleID}_Alt-3-inner.bed"
    ID_Data_3[["chr","start_5", "end_5", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_3_out1,sep="\t",header=False,index=False)
    ID_Data_3[["chr","start_3", "end_3", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_3_out2,sep="\t",header=False,index=False)
    ID_Data_3[["chr","start", "end", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_3_out3,sep="\t",header=False,index=False)
    #############################################################
    alt3_Script1 = f"bedtools intersect -wao -a {ID_Data_3_out1} -b {bedfiledir}/exon-3-start.bed > {outdir}/{SampleID}_Alt-3-5pos-intersect-3.bed"
    alt3_Script2 = f"bedtools intersect -wao -a {ID_Data_3_out2} -b {bedfiledir}/exon-5-start.bed > {outdir}/{SampleID}_Alt-3-3pos-intersect-5.bed"
    alt3_Script3 = f"bedtools intersect -wao -a {ID_Data_3_out3} -b {bedfiledir}/intron.bed > {outdir}/{SampleID}_Alt-3-inner-intersect-intron.bed"
    alt3_Script4 = f"bedtools intersect -wao -a {ID_Data_3_out3} -b {bedfiledir}/exon.bed > {outdir}/{SampleID}_Alt-3-inner-intersect-exon.bed"
    os.system(alt3_Script1)
    os.system(alt3_Script2)
    os.system(alt3_Script3)
    os.system(alt3_Script4)
    Alt_3_Filter_Res = Alt_3_bed_Filter(f"{outdir}/{SampleID}_Alt-3-5pos-intersect-3.bed",
                                        f"{outdir}/{SampleID}_Alt-3-3pos-intersect-5.bed",
                                        f"{outdir}/{SampleID}_Alt-3-inner-intersect-intron.bed",
                                        f"{outdir}/{SampleID}_Alt-3-inner-intersect-exon.bed")
    return Alt_3_Filter_Res

def Alt_3_bed_Filter(alt_splice_unchange,alt_splice_change,alt_splice_inner_intron,alt_splice_inner_exon):
    alt_splice_unchange_data = pd.read_csv(alt_splice_unchange, sep="\t", header=None, skiprows=0)
    alt_splice_change_data = pd.read_csv(alt_splice_change, sep="\t", header=None, skiprows=0)
    alt_splice_inner_intron_data = pd.read_csv(alt_splice_inner_intron, sep="\t", header=None, skiprows=0)
    alt_splice_inner_exon_data = pd.read_csv(alt_splice_inner_exon, sep="\t", header=None, skiprows=0)
    ########################################
    ########################################
    #print(alt_splice_unchange_data)
    alt_splice_unchange_data["gene_Examined"] = alt_splice_unchange_data[3].str.split(":",expand=True)[0]
    alt_splice_unchange_data[["gene_anno","trans_anno"]] = alt_splice_unchange_data[8].str.split(":",expand=True)
    alt_splice_unchange_data_samegene = alt_splice_unchange_data[alt_splice_unchange_data["gene_Examined"]==alt_splice_unchange_data["gene_anno"]]
    alt_splice_unchange_data_samegene = alt_splice_unchange_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #alt_splice_unchange_data_samegene = alt_splice_unchange_data_samegene.groupby([3,9])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    ########################################
    ########################################
    #print(alt_splice_change_data)
    alt_splice_change_data = alt_splice_change_data[[3,10]]
    alt_splice_change_data = alt_splice_change_data.drop_duplicates()
    #alt_splice_change_data.columns = ["Alt_Splice","ChangeSite"]
    #print(alt_splice_change_data)
    ########################################
    ########################################
    #print(alt_splice_inner_intron_data)
    alt_splice_inner_intron_data["len"] = alt_splice_inner_intron_data[7] - alt_splice_inner_intron_data[6]
    alt_splice_inner_intron_data["intersect_ratio"] = alt_splice_inner_intron_data[10] / alt_splice_inner_intron_data["len"]
    alt_splice_inner_intron_data["gene_Examined"] = alt_splice_inner_intron_data[3].str.split(":",expand=True)[0]
    alt_splice_inner_intron_data[["gene_anno","trans_anno"]] = alt_splice_inner_intron_data[8].str.split(":",expand=True)
    alt_splice_inner_intron_data_samegene = alt_splice_inner_intron_data[alt_splice_inner_intron_data["gene_Examined"]==alt_splice_inner_intron_data["gene_anno"]]
    alt_splice_inner_intron_data_samegene = alt_splice_inner_intron_data_samegene[alt_splice_inner_intron_data_samegene["intersect_ratio"]==1]
    alt_splice_inner_intron_data_samegene = alt_splice_inner_intron_data_samegene.groupby([3,"intersect_ratio"])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    ########################################
    ########################################
    #print(alt_splice_inner_exon_data)
    alt_splice_inner_exon_data["len"] = alt_splice_inner_exon_data[7] - alt_splice_inner_exon_data[6]
    alt_splice_inner_exon_data["intersect_ratio"] = alt_splice_inner_exon_data[10] / alt_splice_inner_exon_data["len"]
    ################update
    alt_splice_inner_exon_data_contain = alt_splice_inner_exon_data[alt_splice_inner_exon_data["intersect_ratio"]==1]
    alt_splice_inner_exon_data = alt_splice_inner_exon_data[alt_splice_inner_exon_data["intersect_ratio"]<1]
    alt_splice_inner_exon_data = alt_splice_inner_exon_data[~alt_splice_inner_exon_data[8].isin(alt_splice_inner_exon_data_contain[8])]
    ################update
    alt_splice_inner_exon_data = alt_splice_inner_exon_data[alt_splice_inner_exon_data[10]>=2]
    alt_splice_inner_exon_data["gene_Examined"] = alt_splice_inner_exon_data[3].str.split(":",expand=True)[0]
    alt_splice_inner_exon_data[["gene_anno","trans_anno"]] = alt_splice_inner_exon_data[8].str.split(":",expand=True)
    alt_splice_inner_exon_data_samegene = alt_splice_inner_exon_data[alt_splice_inner_exon_data["gene_Examined"]==alt_splice_inner_exon_data["gene_anno"]]
    alt_splice_inner_exon_data_samegene = alt_splice_inner_exon_data_samegene.groupby([3])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    alt_splice_inner_exon_data_samegene["ExonIntersect"] = 2
    ########################################
    ########################################
    alt_splice_change_data.columns = ["Alt_Splice","ChangeInfo"]
    alt_splice_unchange_data_samegene.columns = ["Alt_Splice","UnChangeInfo","UnChangeTrans"]
    alt_splice_inner_intron_data_samegene.columns = ["Alt_Splice","IntronIntersect","IntronTrans"]
    alt_splice_inner_exon_data_samegene.columns = ["Alt_Splice","ExonTrans","ExonIntersect"]
    alt_splice_change_data = pd.merge(alt_splice_change_data,alt_splice_unchange_data_samegene,on='Alt_Splice',how='left')
    alt_splice_change_data = pd.merge(alt_splice_change_data,alt_splice_inner_intron_data_samegene,on='Alt_Splice',how='left')
    alt_splice_change_data = pd.merge(alt_splice_change_data,alt_splice_inner_exon_data_samegene,on='Alt_Splice',how='left')
    alt_splice_change_data = alt_splice_change_data.fillna(0)
    alt_splice_change_data["Report"] = "YES"
    alt_splice_change_data["TranscriptDetail"] = "YES"
    #outdataframe = outdataframe.head(10)
    for i in alt_splice_change_data.index:
        UnChangeTrans = re.split(r"\|",str(alt_splice_change_data.loc[i,"UnChangeTrans"]))
        IntronTrans = re.split(r"\|",str(alt_splice_change_data.loc[i,"IntronTrans"]))
        ExonTrans = re.split(r"\|",str(alt_splice_change_data.loc[i,"ExonTrans"]))
        UnChangeInfo = alt_splice_change_data.loc[i,"UnChangeInfo"]
        ChangeInfo = alt_splice_change_data.loc[i,"ChangeInfo"]
        common_elements = set(UnChangeTrans) & set(IntronTrans) & set(ExonTrans)
        common_elements = list(common_elements)
        if "0" in common_elements:
            common_elements.remove("0") 
        if len(common_elements) > 0:
            if UnChangeInfo == 1:
                if ChangeInfo == 0:
                    Trans = "YES"
                else:
                    Trans = "NO"
            else:
                Trans = "NO"
            #Trans = "YES"
            TransDetail = "|".join(common_elements)
        else:
            Trans = "NO"
            TransDetail = "-"
        alt_splice_change_data.loc[i,"Report"] = Trans
        alt_splice_change_data.loc[i,"TranscriptDetail"] = TransDetail
    return alt_splice_change_data
    ########################################

def Alt_5_bed(alt5_data,outdir,SampleID,bedfiledir):
    ID_Data_5_out1 = f"{outdir}/{SampleID}_Alt-5-5pos.bed"
    ID_Data_5_out2 = f"{outdir}/{SampleID}_Alt-5-3pos.bed"
    ID_Data_5_out3 = f"{outdir}/{SampleID}_Alt-5-inner.bed"
    alt5_data[["chr","start_5", "end_5", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_5_out1,sep="\t",header=False,index=False)
    alt5_data[["chr","start_3", "end_3", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_5_out2,sep="\t",header=False,index=False)
    alt5_data[["chr","start", "end", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_5_out3,sep="\t",header=False,index=False)
    alt5_Script1 = f"bedtools intersect -wao -a {ID_Data_5_out1} -b {bedfiledir}/exon-3-start.bed > {outdir}/{SampleID}_Alt-5-5pos-intersect-3.bed"
    alt5_Script2 = f"bedtools intersect -wao -a {ID_Data_5_out2} -b {bedfiledir}/exon-5-start.bed > {outdir}/{SampleID}_Alt-5-3pos-intersect-5.bed"
    alt5_Script3 = f"bedtools intersect -wao -a {ID_Data_5_out3} -b {bedfiledir}/intron.bed > {outdir}/{SampleID}_Alt-5-inner-intersect-intron.bed"
    alt5_Script4 = f"bedtools intersect -wao -a {ID_Data_5_out3} -b {bedfiledir}/exon.bed > {outdir}/{SampleID}_Alt-5-inner-intersect-exon.bed"
    os.system(alt5_Script1)
    os.system(alt5_Script2)
    os.system(alt5_Script3)
    os.system(alt5_Script4)
    Alt_5_Filter_Res = Alt_5_bed_Filter(f"{outdir}/{SampleID}_Alt-5-3pos-intersect-5.bed",
                                        f"{outdir}/{SampleID}_Alt-5-5pos-intersect-3.bed",
                                        f"{outdir}/{SampleID}_Alt-5-inner-intersect-intron.bed",
                                        f"{outdir}/{SampleID}_Alt-5-inner-intersect-exon.bed")
    return Alt_5_Filter_Res

def Alt_5_bed_Filter(alt_splice_unchange,alt_splice_change,alt_splice_inner_intron,alt_splice_inner_exon):
    alt_splice_unchange_data = pd.read_csv(alt_splice_unchange, sep="\t", header=None, skiprows=0)
    alt_splice_change_data = pd.read_csv(alt_splice_change, sep="\t", header=None, skiprows=0)
    alt_splice_inner_intron_data = pd.read_csv(alt_splice_inner_intron, sep="\t", header=None, skiprows=0)
    alt_splice_inner_exon_data = pd.read_csv(alt_splice_inner_exon, sep="\t", header=None, skiprows=0)
    ########################################
    ########################################
    #print(alt_splice_unchange_data)
    alt_splice_unchange_data["gene_Examined"] = alt_splice_unchange_data[3].str.split(":",expand=True)[0]
    alt_splice_unchange_data[["gene_anno","trans_anno"]] = alt_splice_unchange_data[8].str.split(":",expand=True)
    alt_splice_unchange_data_samegene = alt_splice_unchange_data[alt_splice_unchange_data["gene_Examined"]==alt_splice_unchange_data["gene_anno"]]
    alt_splice_unchange_data_samegene = alt_splice_unchange_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    alt_splice_unchange_data_samegene = alt_splice_unchange_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    ########################################
    ########################################
    #print(alt_splice_change_data)
    alt_splice_change_data = alt_splice_change_data[[3,10]]
    alt_splice_change_data = alt_splice_change_data.drop_duplicates()
    #alt_splice_change_data.columns = ["Alt_Splice","ChangeSite"]
    #print(alt_splice_change_data)
    ########################################
    ########################################
    #print(alt_splice_inner_intron_data)
    alt_splice_inner_intron_data["len"] = alt_splice_inner_intron_data[7] - alt_splice_inner_intron_data[6]
    alt_splice_inner_intron_data["intersect_ratio"] = alt_splice_inner_intron_data[10] / alt_splice_inner_intron_data["len"]
    alt_splice_inner_intron_data["gene_Examined"] = alt_splice_inner_intron_data[3].str.split(":",expand=True)[0]
    alt_splice_inner_intron_data[["gene_anno","trans_anno"]] = alt_splice_inner_intron_data[8].str.split(":",expand=True)
    alt_splice_inner_intron_data_samegene = alt_splice_inner_intron_data[alt_splice_inner_intron_data["gene_Examined"]==alt_splice_inner_intron_data["gene_anno"]]
    alt_splice_inner_intron_data_samegene = alt_splice_inner_intron_data_samegene[alt_splice_inner_intron_data_samegene["intersect_ratio"]==1]
    alt_splice_inner_intron_data_samegene = alt_splice_inner_intron_data_samegene.groupby([3,"intersect_ratio"])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    ########################################
    ########################################
    #print(alt_splice_inner_exon_data)
    alt_splice_inner_exon_data["len"] = alt_splice_inner_exon_data[7] - alt_splice_inner_exon_data[6]
    alt_splice_inner_exon_data["intersect_ratio"] = alt_splice_inner_exon_data[10] / alt_splice_inner_exon_data["len"]
    ################update
    alt_splice_inner_exon_data_contain = alt_splice_inner_exon_data[alt_splice_inner_exon_data["intersect_ratio"]==1]
    alt_splice_inner_exon_data = alt_splice_inner_exon_data[alt_splice_inner_exon_data["intersect_ratio"]<1]
    alt_splice_inner_exon_data = alt_splice_inner_exon_data[~alt_splice_inner_exon_data[8].isin(alt_splice_inner_exon_data_contain[8])]
    ################update
    alt_splice_inner_exon_data = alt_splice_inner_exon_data[alt_splice_inner_exon_data[10]>=2]
    alt_splice_inner_exon_data["gene_Examined"] = alt_splice_inner_exon_data[3].str.split(":",expand=True)[0]
    alt_splice_inner_exon_data[["gene_anno","trans_anno"]] = alt_splice_inner_exon_data[8].str.split(":",expand=True)
    alt_splice_inner_exon_data_samegene = alt_splice_inner_exon_data[alt_splice_inner_exon_data["gene_Examined"]==alt_splice_inner_exon_data["gene_anno"]]
    alt_splice_inner_exon_data_samegene = alt_splice_inner_exon_data_samegene.groupby([3])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    alt_splice_inner_exon_data_samegene["ExonIntersect"] = 2
    ########################################
    ########################################
    alt_splice_change_data.columns = ["Alt_Splice","ChangeInfo"]
    alt_splice_unchange_data_samegene.columns = ["Alt_Splice","UnChangeInfo","UnChangeTrans"]
    alt_splice_inner_intron_data_samegene.columns = ["Alt_Splice","IntronIntersect","IntronTrans"]
    alt_splice_inner_exon_data_samegene.columns = ["Alt_Splice","ExonTrans","ExonIntersect"]
    alt_splice_change_data = pd.merge(alt_splice_change_data,alt_splice_unchange_data_samegene,on='Alt_Splice',how='left')
    alt_splice_change_data = pd.merge(alt_splice_change_data,alt_splice_inner_intron_data_samegene,on='Alt_Splice',how='left')
    alt_splice_change_data = pd.merge(alt_splice_change_data,alt_splice_inner_exon_data_samegene,on='Alt_Splice',how='left')
    alt_splice_change_data = alt_splice_change_data.fillna(0)
    alt_splice_change_data["Report"] = "YES"
    alt_splice_change_data["TranscriptDetail"] = "YES"
    #outdataframe = outdataframe.head(10)
    for i in alt_splice_change_data.index:
        UnChangeTrans = re.split(r"\|",str(alt_splice_change_data.loc[i,"UnChangeTrans"]))
        IntronTrans = re.split(r"\|",str(alt_splice_change_data.loc[i,"IntronTrans"]))
        ExonTrans = re.split(r"\|",str(alt_splice_change_data.loc[i,"ExonTrans"]))
        UnChangeInfo = alt_splice_change_data.loc[i,"UnChangeInfo"]
        ChangeInfo = alt_splice_change_data.loc[i,"ChangeInfo"]
        common_elements = set(UnChangeTrans) & set(IntronTrans) & set(ExonTrans)
        common_elements = list(common_elements)
        if "0" in common_elements:
            common_elements.remove("0") 
        if len(common_elements) > 0:
            if UnChangeInfo == 1:
                if ChangeInfo == 0:
                    Trans = "YES"
                else:
                    Trans = "NO"
            else:
                Trans = "NO"
            TransDetail = "|".join(common_elements)
        else:
            Trans = "NO"
            TransDetail = "-"
        alt_splice_change_data.loc[i,"Report"] = Trans
        alt_splice_change_data.loc[i,"TranscriptDetail"] = TransDetail
    ########################################
    return alt_splice_change_data
    ########################################

def Alt_C_term_bed(altC_data,outdir,SampleID,bedfiledir):
    ID_Data_altC_out1 = f"{outdir}/{SampleID}_Alt-C-5pos.bed"
    ID_Data_altC_out2 = f"{outdir}/{SampleID}_Alt-C-3pos.bed"
    ID_Data_altC_out3 = f"{outdir}/{SampleID}_Alt-C-inner.bed"
    altC_data[["chr","start_5", "end_5", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altC_out1,sep="\t",header=False,index=False)
    altC_data[["chr","start_3", "end_3", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altC_out2,sep="\t",header=False,index=False)
    altC_data[["chr","start", "end", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altC_out3,sep="\t",header=False,index=False)
    altC_Script1 = f"bedtools intersect -wao -a {ID_Data_altC_out1} -b {bedfiledir}/exon-3-start.bed > {outdir}/{SampleID}_Alt-C-5pos-intersect-3.bed"
    altC_Script2 = f"bedtools intersect -wao -a {ID_Data_altC_out2} -b {bedfiledir}/exon-5-start.bed > {outdir}/{SampleID}_Alt-C-3pos-intersect-5.bed"
    altC_Script3 = f"bedtools intersect -wao -a {ID_Data_altC_out3} -b {bedfiledir}/exon_transcript_last-exon.bed > {outdir}/{SampleID}_Alt-C-inner-intersect-lastExon.bed"
    altC_Script4 = f"bedtools intersect -wao -a {ID_Data_altC_out3} -b {bedfiledir}/exon.bed > {outdir}/{SampleID}_Alt-C-inner-intersect-exon.bed"
    os.system(altC_Script1)
    os.system(altC_Script2)
    os.system(altC_Script3)
    os.system(altC_Script4)
    Alt_C_term_Filter_Res = Alt_C_term_filter(f"{outdir}/{SampleID}_Alt-C-5pos-intersect-3.bed",
                                        f"{outdir}/{SampleID}_Alt-C-3pos-intersect-5.bed",
                                        f"{outdir}/{SampleID}_Alt-C-inner-intersect-lastExon.bed",
                                        f"{outdir}/{SampleID}_Alt-C-inner-intersect-exon.bed",
                                        f"{bedfiledir}/exon_gene_last-exon.bed")
    return Alt_C_term_Filter_Res

def Alt_C_term_filter(alt_splice_unchange,alt_splice_change,alt_splice_last_exon,alt_splice_exon,gene_last_exon):
    alt_splice_unchange_data = pd.read_csv(alt_splice_unchange, sep="\t", header=None, skiprows=0)
    alt_splice_change_data = pd.read_csv(alt_splice_change, sep="\t", header=None, skiprows=0)
    alt_splice_last_exon_data = pd.read_csv(alt_splice_last_exon, sep="\t", header=None, skiprows=0)
    alt_splice_exon_data = pd.read_csv(alt_splice_exon, sep="\t", header=None, skiprows=0)
    gene_last_exon_data = pd.read_csv(gene_last_exon, sep="\t", header=None, skiprows=0)
    #########################################################################################
    #outdataframe = alt_splice_unchange_data[[3]]
    #outdataframe = outdataframe.drop_duplicates()
    outdataframe = alt_splice_unchange_data[[3,4]]
    outdataframe = outdataframe.drop_duplicates()
    outdataframe.columns = ["Alt_Splice","pos_igv_bak"]
    #########################################################################################
    #########################################################################################
    #print(alt_splice_unchange_data)
    alt_splice_unchange_data["gene_Examined"] = alt_splice_unchange_data[3].str.split(":",expand=True)[0]
    alt_splice_unchange_data[["gene_anno","trans_anno"]] = alt_splice_unchange_data[8].str.split(":",expand=True)
    alt_splice_unchange_data_samegene = alt_splice_unchange_data[alt_splice_unchange_data["gene_Examined"]==alt_splice_unchange_data["gene_anno"]]
    alt_splice_unchange_data_samegene = alt_splice_unchange_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #print(alt_splice_change_data)
    #########################################################################################
    #print(alt_splice_unchange_data_samegene)
    #print(alt_splice_change_data_samegene)
    gene_last_exon_data[["gene_anno","trans"]] = gene_last_exon_data[3].str.split(":",expand=True)
    gene_last_exon_data = gene_last_exon_data[["gene_anno",0,1,2,4]]
    gene_last_exon_data.columns = ["gene_anno","chr","start","end","strand"]
    #########################################################################################
    alt_splice_change_data["gene_Examined"] = alt_splice_change_data[3].str.split(":",expand=True)[0]
    alt_splice_change_data[["gene_anno","trans_anno"]] = alt_splice_change_data[8].str.split(":",expand=True)
    alt_splice_change_data_samegene = alt_splice_change_data[alt_splice_change_data["gene_Examined"]==alt_splice_change_data["gene_anno"]]
    alt_splice_change_data_samegene = alt_splice_change_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #########################################################################################
    #print(alt_splice_unchange_data_samegene)
    #print(alt_splice_change_data_samegene)
    #########################################################################################
    #print(alt_splice_last_exon_data)
    alt_splice_last_exon_data["len"] = alt_splice_last_exon_data[7] - alt_splice_last_exon_data[6]
    alt_splice_last_exon_data["intersect_ratio"] = alt_splice_last_exon_data[10] / alt_splice_last_exon_data["len"]
    alt_splice_last_exon_data["gene_Examined"] = alt_splice_last_exon_data[3].str.split(":",expand=True)[0]
    alt_splice_last_exon_data[["gene_anno","trans_anno"]] = alt_splice_last_exon_data[8].str.split(":",expand=True)
    alt_splice_last_exon_data_samegene = alt_splice_last_exon_data[alt_splice_last_exon_data["gene_Examined"]==alt_splice_last_exon_data["gene_anno"]]
    alt_splice_last_exon_data_samegene = alt_splice_last_exon_data_samegene[alt_splice_last_exon_data_samegene["intersect_ratio"]==1]
    alt_splice_last_exon_data_samegene = alt_splice_last_exon_data_samegene.groupby([3,"intersect_ratio"])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #########################################################################################
    #print(alt_splice_exon_data)
    alt_splice_exon_data["gene_Examined"] = alt_splice_exon_data[3].str.split(":",expand=True)[0]
    alt_splice_exon_data[["gene_anno","trans_anno"]] = alt_splice_exon_data[8].str.split(":",expand=True)
    alt_splice_exon_data = alt_splice_exon_data[alt_splice_exon_data["gene_Examined"]==alt_splice_exon_data["gene_anno"]]
    alt_splice_exon_data = alt_splice_exon_data[alt_splice_exon_data[10]>1]
    alt_splice_exon_data = alt_splice_exon_data.groupby([3])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #print(alt_splice_exon_data)
    #########################################################################################
    #print(outdataframe)
    #outdataframe.columns = ["Alt_Splice"]
    alt_splice_unchange_data_samegene.columns = ["Alt_Splice","UnChangeInfo","UnChangeTrans"]
    alt_splice_change_data_samegene.columns = ["Alt_Splice","ChangeInfo","ChangeTrans"]
    alt_splice_last_exon_data_samegene.columns = ["Alt_Splice","LastExonInfo","LastExonTrans"]
    alt_splice_exon_data.columns = ["Alt_Splice","ExonTrans"]
    outdataframe = pd.merge(outdataframe,alt_splice_unchange_data_samegene,on="Alt_Splice",how='left')
    outdataframe = pd.merge(outdataframe,alt_splice_change_data_samegene,on="Alt_Splice",how='left')
    outdataframe = pd.merge(outdataframe,alt_splice_last_exon_data_samegene,on="Alt_Splice",how='left')
    outdataframe = pd.merge(outdataframe,alt_splice_exon_data,on="Alt_Splice",how='left')
    outdataframe["gene_anno"] = outdataframe["Alt_Splice"].str.split(":",expand=True)[0]
    outdataframe = pd.merge(outdataframe,gene_last_exon_data,on="gene_anno",how='left')
    outdataframe = outdataframe.fillna(0)
    outdataframe[["chr_exam","start_exam","end_exam"]] = outdataframe["pos_igv_bak"].str.split(":|-",expand=True)
    #print(outdataframe)
    #########################################################################################
    outdataframe["Report"] = "NO"
    outdataframe["TranscriptDetail"] = "-"
    outdataframe["lastExonDistance"] = "-"
    for i in outdataframe.index:
        Trans = "NO"
        TransDetail = "-"
        UnChangeTrans = re.split(r"\|",str(outdataframe.loc[i,"UnChangeTrans"]))
        UnChangeInfo = outdataframe.loc[i,"UnChangeInfo"]
        ChangeTrans = re.split(r"\|",str(outdataframe.loc[i,"ChangeTrans"]))
        ChangeInfo = outdataframe.loc[i,"ChangeInfo"]
        LastExonTrans = re.split(r"\|",str(outdataframe.loc[i,"LastExonTrans"]))
        LastExonInfo = outdataframe.loc[i,"LastExonInfo"]
        ExonTrans = re.split(r"\|",str(outdataframe.loc[i,"ExonTrans"]))
        ##########################################################
        strand = outdataframe.loc[i,"strand"]
        start = outdataframe.loc[i,"start"]
        end = outdataframe.loc[i,"end"]
        start_exam = outdataframe.loc[i,"start_exam"]
        end_exam = outdataframe.loc[i,"end_exam"]
        if strand == "-":
            distance = int(end) - int(start_exam)
        else:
            distance = int(end_exam) - int(start)
        ##########################################################
        if UnChangeInfo == 1:
            if ChangeInfo == 0:
                if LastExonInfo == 1:
                    common_elements = set(UnChangeTrans) & set(LastExonTrans)
                    if len(common_elements) > 0:
                        diff_elements = set(ChangeTrans) - common_elements
                        if len(diff_elements) > 0:
                            Trans = "YES"
                            TransDetail = "|".join(common_elements) + ":" + "|".join(diff_elements)
                    else:
                        if distance > 0:
                            Trans = "YES"
                            TransDetail = "3UTR"
                if LastExonInfo == 0:
                    common_elements = set(UnChangeTrans) & set(ExonTrans) & set(ChangeTrans)
                    if len(common_elements) > 0:
                        Trans = "YES"
                        TransDetail = "|".join(common_elements)
                    else:
                        if distance > 0:
                            Trans = "YES"
                            TransDetail = "3UTR"
        outdataframe.loc[i,"Report"] = Trans
        outdataframe.loc[i,"TranscriptDetail"] = TransDetail
        outdataframe.loc[i,"lastExonDistance"] = distance
    return outdataframe

def Alt_Promoter_bed(altP_data,outdir,SampleID,bedfiledir):
    ID_Data_altP_out1 = f"{outdir}/{SampleID}_Alt-P-5pos.bed"
    ID_Data_altP_out2 = f"{outdir}/{SampleID}_Alt-P-3pos.bed"
    ID_Data_altP_out3 = f"{outdir}/{SampleID}_Alt-P-inner.bed"
    altP_data[["chr","start_5", "end_5", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altP_out1,sep="\t",header=False,index=False)
    altP_data[["chr","start_3", "end_3", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altP_out2,sep="\t",header=False,index=False)
    altP_data[["chr","start", "end", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altP_out3,sep="\t",header=False,index=False)
    altP_Script1 = f"bedtools intersect -wao -a {ID_Data_altP_out1} -b {bedfiledir}/exon-3-start.bed > {outdir}/{SampleID}_Alt-Promoter-5pos-intersect-3.bed"
    altP_Script2 = f"bedtools intersect -wao -a {ID_Data_altP_out2} -b {bedfiledir}/exon-5-start.bed > {outdir}/{SampleID}_Alt-Promoter-3pos-intersect-5.bed"
    altP_Script3 = f"bedtools intersect -wao -a {ID_Data_altP_out3} -b {bedfiledir}/exon_transcript_first-exon.bed > {outdir}/{SampleID}_Alt-Promoter-inner-intersect-firstexon.bed"
    altP_Script4 = f"bedtools intersect -wao -a {ID_Data_altP_out3} -b {bedfiledir}/exon.bed > {outdir}/{SampleID}_Alt-Promoter-inner-intersect-exon.bed"
    os.system(altP_Script1)
    os.system(altP_Script2)
    os.system(altP_Script3)
    os.system(altP_Script4)
    Alt_Promoter_Filter_Res = Alt_Promoter_bed_Filter(f"{outdir}/{SampleID}_Alt-Promoter-3pos-intersect-5.bed",
                                        f"{outdir}/{SampleID}_Alt-Promoter-5pos-intersect-3.bed",
                                        f"{outdir}/{SampleID}_Alt-Promoter-inner-intersect-firstexon.bed",
                                        f"{outdir}/{SampleID}_Alt-Promoter-inner-intersect-exon.bed",
                                        f"{bedfiledir}/exon_gene_first-exon.bed")
    return Alt_Promoter_Filter_Res

def Alt_Promoter_bed_Filter(alt_splice_unchange,alt_splice_change,alt_splice_first_exon,alt_splice_exon,gene_first_exon):
    alt_splice_unchange_data = pd.read_csv(alt_splice_unchange, sep="\t", header=None, skiprows=0)
    alt_splice_change_data = pd.read_csv(alt_splice_change, sep="\t", header=None, skiprows=0)
    alt_splice_first_exon_data = pd.read_csv(alt_splice_first_exon, sep="\t", header=None, skiprows=0)
    alt_splice_exon_data = pd.read_csv(alt_splice_exon, sep="\t", header=None, skiprows=0)
    gene_first_exon_data = pd.read_csv(gene_first_exon, sep="\t", header=None, skiprows=0)
    #########################################################################################
    #outdataframe = alt_splice_unchange_data[[3]]
    #outdataframe = outdataframe.drop_duplicates()
    outdataframe = alt_splice_unchange_data[[3,4]]
    outdataframe = outdataframe.drop_duplicates()
    outdataframe.columns = ["Alt_Splice","pos_igv_bak"]
    #########################################################################################
    #print(alt_splice_unchange_data)
    alt_splice_unchange_data["gene_Examined"] = alt_splice_unchange_data[3].str.split(":",expand=True)[0]
    alt_splice_unchange_data[["gene_anno","trans_anno"]] = alt_splice_unchange_data[8].str.split(":",expand=True)
    alt_splice_unchange_data_samegene = alt_splice_unchange_data[alt_splice_unchange_data["gene_Examined"]==alt_splice_unchange_data["gene_anno"]]
    alt_splice_unchange_data_samegene = alt_splice_unchange_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #print(alt_splice_change_data)
    #########################################################################################
    alt_splice_change_data["gene_Examined"] = alt_splice_change_data[3].str.split(":",expand=True)[0]
    alt_splice_change_data[["gene_anno","trans_anno"]] = alt_splice_change_data[8].str.split(":",expand=True)
    alt_splice_change_data_samegene = alt_splice_change_data[alt_splice_change_data["gene_Examined"]==alt_splice_change_data["gene_anno"]]
    alt_splice_change_data_samegene = alt_splice_change_data_samegene.groupby([3,10])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #########################################################################################
    #print(alt_splice_unchange_data_samegene)
    #print(alt_splice_change_data_samegene)
    gene_first_exon_data[["gene_anno","trans"]] = gene_first_exon_data[3].str.split(":",expand=True)
    gene_first_exon_data = gene_first_exon_data[["gene_anno",0,1,2,4]]
    gene_first_exon_data.columns = ["gene_anno","chr","start","end","strand"]
    #########################################################################################
    #print(alt_splice_last_exon_data)
    alt_splice_first_exon_data["len"] = alt_splice_first_exon_data[7] - alt_splice_first_exon_data[6]
    alt_splice_first_exon_data["intersect_ratio"] = alt_splice_first_exon_data[10] / alt_splice_first_exon_data["len"]
    alt_splice_first_exon_data["gene_Examined"] = alt_splice_first_exon_data[3].str.split(":",expand=True)[0]
    alt_splice_first_exon_data[["gene_anno","trans_anno"]] = alt_splice_first_exon_data[8].str.split(":",expand=True)
    alt_splice_first_exon_data_samegene = alt_splice_first_exon_data[alt_splice_first_exon_data["gene_Examined"]==alt_splice_first_exon_data["gene_anno"]]
    alt_splice_first_exon_data_samegene = alt_splice_first_exon_data_samegene[alt_splice_first_exon_data_samegene["intersect_ratio"]==1]
    alt_splice_first_exon_data_samegene = alt_splice_first_exon_data_samegene.groupby([3,"intersect_ratio"])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #########################################################################################
    #print(alt_splice_exon_data)
    alt_splice_exon_data["gene_Examined"] = alt_splice_exon_data[3].str.split(":",expand=True)[0]
    alt_splice_exon_data[["gene_anno","trans_anno"]] = alt_splice_exon_data[8].str.split(":",expand=True)
    alt_splice_exon_data = alt_splice_exon_data[alt_splice_exon_data["gene_Examined"]==alt_splice_exon_data["gene_anno"]]
    alt_splice_exon_data = alt_splice_exon_data[alt_splice_exon_data[10]>1]
    alt_splice_exon_data = alt_splice_exon_data.groupby([3])["trans_anno"].agg(lambda x: '|'.join(x)).reset_index()
    #########################################################################################
    #print(outdataframe)
    #outdataframe.columns = ["Alt_Splice"]
    alt_splice_unchange_data_samegene.columns = ["Alt_Splice","UnChangeInfo","UnChangeTrans"]
    alt_splice_change_data_samegene.columns = ["Alt_Splice","ChangeInfo","ChangeTrans"]
    alt_splice_first_exon_data_samegene.columns = ["Alt_Splice","FirstExonInfo","FirstExonTrans"]
    alt_splice_exon_data.columns = ["Alt_Splice","ExonTrans"]
    outdataframe = pd.merge(outdataframe,alt_splice_unchange_data_samegene,on="Alt_Splice",how='left')
    outdataframe = pd.merge(outdataframe,alt_splice_change_data_samegene,on="Alt_Splice",how='left')
    outdataframe = pd.merge(outdataframe,alt_splice_first_exon_data_samegene,on="Alt_Splice",how='left')
    outdataframe = pd.merge(outdataframe,alt_splice_exon_data,on="Alt_Splice",how='left')
    outdataframe["gene_anno"] = outdataframe["Alt_Splice"].str.split(":",expand=True)[0]
    outdataframe = pd.merge(outdataframe,gene_first_exon_data,on="gene_anno",how='left')
    outdataframe = outdataframe.fillna(0)
    outdataframe[["chr_exam","start_exam","end_exam"]] = outdataframe["pos_igv_bak"].str.split(":|-",expand=True)
    print(outdataframe)
    #########################################################################################
    outdataframe["Report"] = "NO"
    outdataframe["TranscriptDetail"] = "-"
    outdataframe["FirstExonDistance"] = "-"
    for i in outdataframe.index:
        Trans = "NO"
        TransDetail = "-"
        UnChangeTrans = re.split(r"\|",str(outdataframe.loc[i,"UnChangeTrans"]))
        UnChangeInfo = outdataframe.loc[i,"UnChangeInfo"]
        ChangeTrans = re.split(r"\|",str(outdataframe.loc[i,"ChangeTrans"]))
        ChangeInfo = outdataframe.loc[i,"ChangeInfo"]
        FirstExonTrans = re.split(r"\|",str(outdataframe.loc[i,"FirstExonTrans"]))
        FirstExonInfo = outdataframe.loc[i,"FirstExonInfo"]
        ExonTrans = re.split(r"\|",str(outdataframe.loc[i,"ExonTrans"]))
        strand = outdataframe.loc[i,"strand"]
        start = outdataframe.loc[i,"start"]
        end = outdataframe.loc[i,"end"]
        start_exam = outdataframe.loc[i,"start_exam"]
        end_exam = outdataframe.loc[i,"end_exam"]
        if strand == "+":
            distance = int(end) - int(start_exam)
        else:
            distance = int(end_exam) - int(start)
        if UnChangeInfo == 1:
            if ChangeInfo == 0:
                if FirstExonInfo == 1:
                    common_elements = set(UnChangeTrans) & set(FirstExonTrans)
                    if len(common_elements) > 0:
                        diff_elements = set(ChangeTrans) - common_elements
                        if len(diff_elements) > 0:
                            Trans = "YES"
                            TransDetail = "|".join(common_elements) + ":" + "|".join(diff_elements)
                    else:
                        if distance > 0:
                            Trans = "YES"
                            TransDetail = "5UTR"
                if FirstExonInfo == 0:
                    common_elements = set(UnChangeTrans) & set(ExonTrans) & set(ChangeTrans)
                    if len(common_elements) > 0:
                        Trans = "YES"
                        TransDetail = "|".join(common_elements)
                    else:
                        if distance > 0:
                            Trans = "YES"
                            TransDetail = "5UTR"
        outdataframe.loc[i,"Report"] = Trans
        outdataframe.loc[i,"TranscriptDetail"] = TransDetail
        outdataframe.loc[i,"FirstExonDistance"] = distance
    return outdataframe

def Trans_Splicing_bed(altT_data,outdir,SampleID,bedfiledir):
    ID_Data_altT_out1 = f"{outdir}/{SampleID}_Alt-TransSplice-5pos.bed"
    ID_Data_altT_out2 = f"{outdir}/{SampleID}_Alt-TransSplice-3pos.bed"
    altT_data[["chr","start_5", "end_5", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altT_out1,sep="\t",header=False,index=False)
    altT_data[["chr","start_3", "end_3", "Alt_Splice","pos_IGV"]].to_csv(ID_Data_altT_out2,sep="\t",header=False,index=False)
    altT_Script1 = f"bedtools intersect -wao -a {ID_Data_altT_out1} -b {bedfiledir}/gene.bed > {outdir}/{SampleID}_Alt-TransSplice-5pos-gene.bed"
    altT_Script2 = f"bedtools intersect -wao -a {ID_Data_altT_out2} -b {bedfiledir}/gene.bed > {outdir}/{SampleID}_Alt-TransSplice-3pos-gene.bed"
    os.system(altT_Script1)
    os.system(altT_Script2)
    Trans_Splicing_Filter_data = Trans_Splicing_bed_Filter(f"{outdir}/{SampleID}_Alt-TransSplice-3pos-gene.bed",
                                        f"{outdir}/{SampleID}_Alt-TransSplice-5pos-gene.bed")
    return Trans_Splicing_Filter_data

def Trans_Splicing_bed_Filter(pos3_intersect_file,pos5_intersect_file):
    #########################################################################################
    pos3_intersect_data = pd.read_csv(pos3_intersect_file, sep="\t", header=None, skiprows=0)
    pos5_intersect_data = pd.read_csv(pos5_intersect_file, sep="\t", header=None, skiprows=0)
    #print(pos3_intersect_data)
    #print(pos5_intersect_data)
    #########################################################################################
    #pos3_intersect_data["gene_Examined"] = pos3_intersect_data[3].str.split(":",expand=True)[0]
    pos3_intersect_data = pos3_intersect_data.groupby([3])[8].agg(lambda x: '|'.join(x)).reset_index()
    pos5_intersect_data = pos5_intersect_data.groupby([3])[8].agg(lambda x: '|'.join(x)).reset_index()
    intersect_data_merge = pd.merge(pos3_intersect_data,pos5_intersect_data,on=3,how='left')
    intersect_data_merge.columns = ["Alt_Splice","pos3_gene","pos5_gene"]
    intersect_data_merge["Report"] = "-"
    for i in intersect_data_merge.index:
        Alt_Splice = intersect_data_merge.loc[i,"Alt_Splice"]
        pos3_gene = re.split(r"\|",str(intersect_data_merge.loc[i,"pos3_gene"]))
        pos5_gene = re.split(r"\|",str(intersect_data_merge.loc[i,"pos5_gene"]))
        common_elements = set(pos3_gene) & set(pos5_gene)
        if len(common_elements) <= 0 & Alt_Splice.find("-ENSG")>=0:
            intersect_data_merge.loc[i,"Report"] = "YES"
        else:
            intersect_data_merge.loc[i,"Report"] = "NO"
        if Alt_Splice.find("-ENSG")<0:
            intersect_data_merge.loc[i,"Report"] = "NO"
    return intersect_data_merge


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--prunedcount", help="prunedcount file", required=True)
    parser.add_argument("-i", "--intronJunction", help="intronJunction file", required=True)
    parser.add_argument("-j", "--junctionbed", help="junctionbed file", required=True)
    parser.add_argument("-e", "--EventAnnotation", help="EventAnnotation file", required=True)
    parser.add_argument("-o", "--outputdir", help="output dir", required=True)
    parser.add_argument("-s", "--SampleID", help="SampleID", required=True)
    parser.add_argument("-g", "--Gtf_Bed_File_Dir", help="Gtf_Bed_File_Dir", required=True)
    parser.add_argument("-b", "--outFilterbam", help="outFilterbam file", required=True)
    args = parser.parse_args()
    outputdir = args.outputdir
    outFilterbam = args.outFilterbam
    output_Bed_Dir = f"{outputdir}/Bed_FileDir"
    igv_dir = outFilterbam.replace("/mnt/user/dingyu/01.Shared_Neoantigen/01.SNAF/05.Clinical_Sample_Analysis","http://58.242.248.155:5003")
    os.makedirs(output_Bed_Dir, exist_ok=True)
    ############################################################################################
    ############################################################################################
    ############################################################################################
    purned_counts_data = pd.read_csv(args.prunedcount,sep="\t",header=0,skiprows=0)
    purned_counts_data = purned_counts_data.rename(columns={'Unnamed: 0': 'Alt_Splice'})
    purned_counts_data = purned_counts_data[["Alt_Splice",f"{args.SampleID}_tumor.bed"]]
    purned_counts_data.columns = ["Alt_Splice","PrunedCount"]
    EventAnnotation = pd.read_csv(args.EventAnnotation,sep="\t",header=0)
    EventAnnotation = EventAnnotation[["Symbol","Examined-Junction","Background-Major-Junction","AltExons","Coordinates","EventAnnotation"]]
    EventAnnotation[["pos_Examined","pos_BKG"]] = EventAnnotation["Coordinates"].str.split("|", expand=True, regex=False)
    EventAnnotation = EventAnnotation.drop(["Coordinates"],axis=1)
    EventAnnotation.rename(columns={'Examined-Junction': 'Alt_Splice'}, inplace=True)
    EventAnnotation.rename(columns={'Background-Major-Junction': 'BKG_Splice'}, inplace=True)
    EventAnnotation = EventAnnotation[EventAnnotation["Alt_Splice"].isin(purned_counts_data["Alt_Splice"])]
    EventAnnotation = pd.merge(EventAnnotation,purned_counts_data,how='inner',on='Alt_Splice')
    EventAnnotation = EventAnnotation[~EventAnnotation["pos_Examined"].str.contains('chrY')]
    #print(EventAnnotation)
    #print(purned_counts_data)
    ############################################################################################
    ############################################################################################
    ####内含子保留注释信息提取
    EventAnnotation_IR = EventAnnotation[EventAnnotation["EventAnnotation"] == "intron-retention"]
    #EventAnnotation_IR = EventAnnotation_IR.fillna(0)
    #EventAnnotation_IR = EventAnnotation_IR[EventAnnotation_IR[f"{args.SampleID}.Aligned.sortedByCoord.out.bed"]>0]
    EventAnnotation_IR["MergeKey"] = EventAnnotation_IR["Alt_Splice"].str.split(':',expand=True)[0] + "_" + EventAnnotation_IR["pos_Examined"]
    ############################################################################################
    ############################################################################################
    ####从intronJunction的bed文件中提取内含子保留事件的信息，用于回溯到底哪两个断点为一个完整的内含子保留事件
    try:
        intronJunction = pd.read_csv(args.intronJunction,sep="\t",header=None)
        intronJunction_Process = bed_file_Process(intronJunction)
        intronJunction_Process_merge = pd.merge(EventAnnotation_IR,intronJunction_Process,how='inner',on='MergeKey')
        pair_info_data_count = intronJunction_Process_merge.groupby(['SpliceInfo']).size().reset_index(name='Count')
        pair_info_data_count = pair_info_data_count.sort_values(by='Count', ascending=False)
        pair_info_data_count = pair_info_data_count[pair_info_data_count["Count"]==2]
        intronJunction_Process_merge = intronJunction_Process_merge[intronJunction_Process_merge["SpliceInfo"].isin(pair_info_data_count["SpliceInfo"])]
        intronJunction_Process_merge = intronJunction_Process_merge.sort_values(by="SpliceInfo",ascending=False)
        intronJunction_Process_merge = intronJunction_Process_merge.sort_values(by=["SpliceInfo","pos_Examined"],ascending=[True,True])
        intronJunction_Process_merge.reset_index(drop=True,inplace=True)
        odd_sample = intronJunction_Process_merge[intronJunction_Process_merge.index % 2 == 0]
        even_sample = intronJunction_Process_merge[intronJunction_Process_merge.index % 2 == 1]
        intronEvent_merge_data = pd.merge(odd_sample,even_sample,on="SpliceInfo",how='inner')
        #print(intronEvent_merge_data)
        #print(intronEvent_merge_data.loc[0])
        #print(intronEvent_merge_data.keys())
    except:
        #intronJunction_Process = pd.DataFrame(columns=["MergeKey","SpliceInfo","pos_Splice","Strand","support_reads_intronbed"])
        intronEvent_merge_data = pd.DataFrame(columns=['Symbol_x', 'Alt_Splice_x', 'BKG_Splice_x', 'AltExons_x',\
                                            'EventAnnotation_x', 'pos_Examined_x', 'pos_BKG_x', 'PrunedCount_x',\
                                            'MergeKey_x', 'SpliceInfo', 'pos_Splice_x', 'Strand_x',\
                                            'support_reads_intronbed_x', 'Symbol_y', 'Alt_Splice_y', 'BKG_Splice_y',\
                                            'AltExons_y', 'EventAnnotation_y', 'pos_Examined_y', 'pos_BKG_y',\
                                            'PrunedCount_y', 'MergeKey_y', 'pos_Splice_y', 'Strand_y',\
                                            'support_reads_intronbed_y'])
    if len(intronEvent_merge_data) > 0:
        intronEvent_merge_data[["chr_x","start_x","end_x"]] = intronEvent_merge_data["pos_Examined_x"].str.split(":|-",expand=True)
        intronEvent_merge_data[["chr_y","start_y","end_y"]] = intronEvent_merge_data["pos_Examined_y"].str.split(":|-",expand=True)
    else:
        intronEvent_merge_data[["chr_x","start_x","end_x","chr_y","start_y","end_y"]] = ""
    intronEvent_merge_data["start_x"] = intronEvent_merge_data["start_x"].astype(int)
    intronEvent_merge_data["end_x"] = intronEvent_merge_data["end_x"].astype(int)
    intronEvent_merge_data["start_y"] = intronEvent_merge_data["start_y"].astype(int)
    intronEvent_merge_data["end_y"] = intronEvent_merge_data["end_y"].astype(int)
    output_data_plus = intronEvent_merge_data[intronEvent_merge_data["Strand_x"]=="+"]
    output_data_plus["pos_IGV"] = output_data_plus["chr_x"] + ":" + (output_data_plus[["start_x","start_y"]].min(axis=1)).astype(str) + "-" + (output_data_plus[["end_x","end_y"]].max(axis=1)).astype(str)
    output_data_minus = intronEvent_merge_data[intronEvent_merge_data["Strand_x"]=="-"]
    output_data_minus["pos_IGV"] = output_data_minus["chr_x"] + ":" + (output_data_minus[["end_x","end_y"]].min(axis=1)).astype(str) + "-" + (output_data_minus[["start_x","start_y"]].max(axis=1)).astype(str)
    intronEvent_merge_Res = pd.concat([output_data_plus,output_data_minus])
    intronEvent_merge_Res["SampleID"] = args.SampleID
    intronEvent_merge_Res["bamDir"] = igv_dir
    IR_Filter_Res = IR_bed(intronEvent_merge_Res,output_Bed_Dir,args.SampleID,args.Gtf_Bed_File_Dir)
    if len(IR_Filter_Res) > 0:
        intronEvent_merge_Res = pd.merge(IR_Filter_Res,intronEvent_merge_Res,how='left',on='SpliceInfo')
    else:
        IR_Filter_Res = pd.DataFrame({
            "SpliceInfo":[],
            "IntronGene":[],
            "IntronTrans":[],
            "exonGene":[],
            "exonTrans":[],
            "Pos1Gene":[],
            "Pos1Trans":[],
            "Pos2Gene":[],
            "Pos2Trans":[],
            "intersect_exonmerge":[],
            "Report":[],
            "GeneDetail":[],
            "TranscriptDetail":[],
        })
        #intronEvent_merge_Res = pd.merge(IR_Filter_Res,intronEvent_merge_Res,how='left',on='SpliceInfo')
        intronEvent_merge_Res = intronEvent_merge_Res.drop(columns=["SpliceInfo"])
        intronEvent_merge_Res = pd.concat([IR_Filter_Res,intronEvent_merge_Res],axis=1)
        #print(intronEvent_merge_Res)
    intronEvent_merge_Res.to_csv(f"{outputdir}/{args.SampleID}_intron-retention_Raw.txt",sep="\t",index=False,header=True)
    #intronEvent_merge_Res[intronEvent_merge_Res["Report"]=="YES"].to_csv(f"{outputdir}/{args.SampleID}_intron-retention_Filter.txt",sep="\t",index=False,header=True)
    #print(intronEvent_merge_Res)
    #print(IR_Filter_Res)
    ############################################################################################
    ############################################################################################
    #print(intronJunction_Process_merge)
    #print(intronEvent_merge_Res)
    ############################################################################################
    ############################################################################################
    #junctionbed = junctionbed.head(5000)
    try:
        junctionbed = pd.read_csv(args.junctionbed,sep="\t",header=None,skiprows=1)
        junctionbed_Process = bed_file_Process_2(junctionbed)
        ####其他可变剪接事件注释信息提取
        EventAnnotation_other = EventAnnotation[EventAnnotation["EventAnnotation"] != "intron-retention"]
        EventAnnotation_other_merge = pd.merge(EventAnnotation_other,junctionbed_Process,how='inner',on='pos_Examined')
        EventAnnotation_other_merge[["chr","start_x","end_x"]] = EventAnnotation_other_merge["pos_Examined"].str.split(":|-",expand=True)
        EventAnnotation_other_merge["start_x"] = EventAnnotation_other_merge["start_x"].astype(int)
        EventAnnotation_other_merge["end_x"] = EventAnnotation_other_merge["end_x"].astype(int)
        EventAnnotation_other_merge["start"] = EventAnnotation_other_merge[["start_x","end_x"]].min(axis=1)
        EventAnnotation_other_merge["end"] = EventAnnotation_other_merge[["start_x","end_x"]].max(axis=1)
        EventAnnotation_other_merge["pos_IGV"] = EventAnnotation_other_merge["chr"] + ":" + (EventAnnotation_other_merge["start"]).astype(str) + "-" + (EventAnnotation_other_merge["end"]).astype(str)
        EventOther_merge = EventAnnotation_other_merge[["Alt_Splice", "Symbol" ,"EventAnnotation", "chr","start","end","start_x","end_x","Strand","pos_Examined","pos_IGV","PrunedCount"]]
        EventOther_merge = EventOther_merge.drop_duplicates()
        EventOther_merge["start_5"] = EventOther_merge["start_x"] -1 
        EventOther_merge["end_5"] = EventOther_merge["start_x"]
        EventOther_merge["start_3"] = EventOther_merge["end_x"] -1 
        EventOther_merge["end_3"] = EventOther_merge["end_x"]
        EventOther_merge["start_left"] = (EventOther_merge["start"].astype(int)) - 1
        EventOther_merge["start_right"] = (EventOther_merge["end"].astype(int)) - 1  
        EventOther_merge = EventOther_merge[["chr","start","end","start_5","end_5","start_3","end_3","start_left","start_right","Alt_Splice","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand"]]
        EventOther_merge = EventOther_merge.drop_duplicates()
    except:
        EventOther_merge = pd.DataFrame(columns=["chr","start","end","start_5","end_5","start_3","end_3","start_left","start_right","Alt_Splice","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand"])
    #print(EventOther_merge)
    ############################################################################################
    ############################################################################################
    alt_type_list = ["cassette-exon","alt-3'","alt-5'","alt-C-term","altPromoter","trans-splicing"]
    for alt_type in alt_type_list:
        EventOther_merge_alt_type = EventOther_merge[EventOther_merge["EventAnnotation"] == alt_type]
        EventOther_merge_alt_type["SampleID"] = args.SampleID
        EventOther_merge_alt_type["bamDir"] = igv_dir
        #print(EventOther_merge_alt_type)
        #if len(EventOther_merge_alt_type) > 0:
        if alt_type == "cassette-exon":
            if len(EventOther_merge_alt_type) > 0:
                Exon_Skip_Res = Exon_Skip_bed(EventOther_merge_alt_type,output_Bed_Dir,args.SampleID,args.Gtf_Bed_File_Dir)
                Exon_Skip_Res = pd.merge(Exon_Skip_Res,EventOther_merge_alt_type,how='left',on='Alt_Splice')
            else:
                Exon_Skip_Res = pd.DataFrame(columns=["Alt_Splice","leftpos_samegene","leftpos_samegene_trans","rightpos_samegene","rightpos_samegene_trans","inner_exon_samegene","inner_exon_samegene_trans","inner_intron_samegene_trans","inner_intron_samegene","Report","TranscriptDetail","chr","start","end","start_5","end_5","start_3","end_3","start_left","start_right","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand","SampleID","bamDir"])
            Exon_Skip_Res.to_csv(f"{outputdir}/{args.SampleID}_cassette-exon_Raw.txt",sep="\t",index=False,header=True)
            #Exon_Skip_Res[Exon_Skip_Res["Report"]=="YES"].to_csv(f"{outputdir}/{args.SampleID}_cassette-exon_Filter.txt",sep="\t",index=False,header=True)
            #print(Exon_Skip_Res)
        if alt_type == "alt-3'":
            if len(EventOther_merge_alt_type) > 0:
                Alt_3_Res = Alt_3_bed(EventOther_merge_alt_type,output_Bed_Dir,args.SampleID,args.Gtf_Bed_File_Dir)
                Alt_3_Res = pd.merge(Alt_3_Res,EventOther_merge_alt_type,how='left',on='Alt_Splice')
            else:
                Alt_3_Res = pd.DataFrame(columns=["Alt_Splice","ChangeInfo","UnChangeInfo","UnChangeTrans","IntronIntersect","IntronTrans","ExonTrans","ExonIntersect","Report","TranscriptDetail","chr","start","end","start_5","end_5","start_3","end_3","start_left","start_right","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand","SampleID","bamDir"])
            Alt_3_Res.to_csv(f"{outputdir}/{args.SampleID}_alt-3_Raw.txt",sep="\t",index=False,header=True)
            #Alt_3_Res[Alt_3_Res["Report"]=="YES"].to_csv(f"{outputdir}/{args.SampleID}_alt-3_Filter.txt",sep="\t",index=False,header=True)
            #print(Alt_3_Res)
        if alt_type == "alt-5'":
            if len(EventOther_merge_alt_type) > 0:
                Alt_5_Res = Alt_5_bed(EventOther_merge_alt_type,output_Bed_Dir,args.SampleID,args.Gtf_Bed_File_Dir)
                Alt_5_Res = pd.merge(Alt_5_Res,EventOther_merge_alt_type,how='left',on='Alt_Splice')
            else:
                Alt_5_Res = pd.DataFrame(columns=["Alt_Splice","ChangeInfo","UnChangeInfo","UnChangeTrans","IntronIntersect","IntronTrans","ExonTrans","ExonIntersect","Report","TranscriptDetail","chr","start","end","start_5","end_5","start_3","end_3","start_left","start_right","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand","SampleID","bamDir"])
            Alt_5_Res.to_csv(f"{outputdir}/{args.SampleID}_alt-5_Raw.txt",sep="\t",index=False,header=True)
            #Alt_5_Res[Alt_5_Res["Report"]=="YES"].to_csv(f"{outputdir}/{args.SampleID}_alt-5_Filter.txt",sep="\t",index=False,header=True)
            #print(Alt_5_Res)
        if alt_type == "alt-C-term":
            if len(EventOther_merge_alt_type) > 0:
                Alt_C_term_Res = Alt_C_term_bed(EventOther_merge_alt_type,output_Bed_Dir,args.SampleID,args.Gtf_Bed_File_Dir)
                Alt_C_term_Res = pd.merge(Alt_C_term_Res,EventOther_merge_alt_type,how='left',on='Alt_Splice')
            else:
                Alt_C_term_Res = pd.DataFrame(columns=["Alt_Splice","pos_igv_bak","UnChangeInfo","UnChangeTrans","ChangeInfo","ChangeTrans","LastExonInfo","LastExonTrans","ExonTrans","gene_anno","chr_x","start_x","end_x","strand","chr_exam","start_exam","end_exam","Report","TranscriptDetail","lastExonDistance","chr_y","start_y","end_y","start_5","end_5","start_3","end_3","start_left","start_right","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand","SampleID","bamDir"])
            Alt_C_term_Res.to_csv(f"{outputdir}/{args.SampleID}_alt-C-term_Raw.txt",sep="\t",index=False,header=True)
            #Alt_C_term_Res[Alt_C_term_Res["Report"]=="YES"].to_csv(f"{outputdir}/{args.SampleID}_alt-C-term_Filter.txt",sep="\t",index=False,header=True)
            #print(Alt_C_term_Res)
        if alt_type == "altPromoter":
            if len(EventOther_merge_alt_type) > 0:
                Alt_Promoter_Res = Alt_Promoter_bed(EventOther_merge_alt_type,output_Bed_Dir,args.SampleID,args.Gtf_Bed_File_Dir)
                Alt_Promoter_Res = pd.merge(Alt_Promoter_Res,EventOther_merge_alt_type,how='left',on='Alt_Splice')
            else:
                Alt_Promoter_Res = pd.DataFrame(columns=["Alt_Splice","pos_igv_bak","UnChangeInfo","UnChangeTrans","ChangeInfo","ChangeTrans","FirstExonInfo","FirstExonTrans","ExonTrans","gene_anno","chr_x","start_x","end_x","strand","chr_exam","start_exam","end_exam","Report","TranscriptDetail","FirstExonDistance","chr_y","start_y","end_y","start_5","end_5","start_3","end_3","start_left","start_right","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand","SampleID","bamDir"])
            Alt_Promoter_Res.to_csv(f"{outputdir}/{args.SampleID}_altPromoter_Raw.txt",sep="\t",index=False,header=True)
            #Alt_Promoter_Res[Alt_Promoter_Res["Report"]=="YES"].to_csv(f"{outputdir}/{args.SampleID}_altPromoter_Filter.txt",sep="\t",index=False,header=True)
            #print(Alt_Promoter_Res)
        if alt_type == "trans-splicing":
            if len(EventOther_merge_alt_type) > 0:
                Trans_Splicing_Res = Trans_Splicing_bed(EventOther_merge_alt_type,output_Bed_Dir,args.SampleID,args.Gtf_Bed_File_Dir)
                Trans_Splicing_Res = pd.merge(Trans_Splicing_Res,EventOther_merge_alt_type,how='left',on='Alt_Splice')
            else:
                Trans_Splicing_Res = pd.DataFrame(columns=["Alt_Splice","pos3_gene","pos5_gene","Report","chr","start","end","start_5","end_5","start_3","end_3","start_left","start_right","pos_Examined","pos_IGV","PrunedCount","EventAnnotation","Symbol","Strand","SampleID","bamDir"])
            Trans_Splicing_Res.to_csv(f"{outputdir}/{args.SampleID}_trans-splicing_Raw.txt",sep="\t",index=False,header=True)
            #Trans_Splicing_Res[Trans_Splicing_Res["Report"]=="YES"].to_csv(f"{outputdir}/{args.SampleID}_trans-splicing_Filter.txt",sep="\t",index=False,header=True)
            #print(Trans_Splicing_Res)

