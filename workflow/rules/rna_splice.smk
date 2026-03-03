# Helper function to get sample names from dataframe
# (You might need to adjust based on how your 'samples' variable is structured)


rule bam_to_bed:
    input:
        bam = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam",
        bai = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam.bai"
    output:
        # We define the expected output files for one sample
        # Note: AltAnalyze script naming convention usually appends '__junction.bed'
        junc_bed = "rna/altanalyze/bed/{sample}_{sampletype}__junction.bed",
        exon_bed = "rna/altanalyze/bed/{sample}_{sampletype}__intronJunction.bed"
    log:
        "logs/altanalyze/bam_to_bed_{sample}_{sampletype}.log"
    params:
        db_path = config["reference"]["AltDatabase"],
        bed_dir = "rna/altanalyze/bed",
        species = "Hs",
        version = "EnsMart91"
    conda:
        "../../envs/altanalyze.yaml"
    shell:
        """
        # 0. SETUP
        TOOL_DIR="$CONDA_PREFIX/share/AltAnalyze"
        
        # Ensure tool is installed (Clone if missing)
        # We need to lock this or assume it's installed to prevent race conditions in parallel
        if [ ! -f "$TOOL_DIR/AltAnalyze.py" ]; then
            # Simple lock mechanism or just fail first time
            echo "Cloning AltAnalyze..." >> {log}
            mkdir -p "$CONDA_PREFIX/share"
            git clone https://github.com/nsalomonis/altanalyze.git "$TOOL_DIR" >> {log} 2>&1
        fi

        # 1. PREPARE INPUTS
        # Create a temp symlink with the exact name we want output to have
        # The script outputs files based on input filename
        
        # We use a unique temp identifier for the link to avoid collisions
        TEMP_BAM="{params.bed_dir}/{wildcards.sample}_{wildcards.sampletype}.bam"
        
        mkdir -p {params.bed_dir}
        ln -sf $(readlink -f {input.bam}) "$TEMP_BAM"
        ln -sf $(readlink -f {input.bai}) "$TEMP_BAM.bai"
        
        echo "Converting $TEMP_BAM..." > {log}

        # 2. RUN CONVERSION
        python "$TOOL_DIR/import_scripts/BAMtoJunctionBED.py" \
            --i "$TEMP_BAM" \
            --species {params.species} \
            --r "{params.db_path}/{params.version}/ensembl/{params.species}/{params.species}_Ensembl_exon.txt" \
            >> {log} 2>&1
            
        python "$TOOL_DIR/import_scripts/BAMtoExonBED.py" \
            --i "$TEMP_BAM" \
            --r "{params.db_path}/{params.version}/ensembl/{params.species}/{params.species}.bed" \
            --s {params.species} \
            >> {log} 2>&1
        """


def get_rna_samples():
    return samples[samples['data_type'] == 'rna']
# Helper to get all expected bed files
def get_all_bed_files(file_type="junction"):
    files = []
    for index, row in get_rna_samples().iterrows():
        s = row["sample_id"]
        t = row["sampletype"]
        files.append(f"rna/altanalyze/bed/{s}_{t}__{file_type}.bed")
    return files

rule AltAnalyze_run:
    input:
        # Require all BED files from the previous bam_to_bed rule
        junc_beds = get_all_bed_files("junction"),
        exon_beds = get_all_bed_files("intronJunction"),
        prune_script = workflow.source_path("../scripts/prune_altanalyze.py")
    output:
        counts_full_file = "rna/altanalyze/output/ExpressionInput/counts.original.full.txt",
        counts_pruned_file = "rna/altanalyze/output/ExpressionInput/counts.original.pruned.txt",
        counts_pruned_spec_file = "rna/altanalyze/output/ExpressionInput/counts.original.specific.txt",

    log:
        "logs/altanalyze/AltAnalyze_main.log"
    params:
        db_path = config["reference"]["AltDatabase"],
        output_dir = "rna/altanalyze/output",
        bed_dir = "rna/altanalyze/bed",
        species = "Hs",
        platform = "RNASeq",
        version = "EnsMart91",
        task = "original",
    conda:
        "../../envs/altanalyze.yaml"
    shell:
        """
        TOOL_DIR="$CONDA_PREFIX/share/AltAnalyze"
        
        # ---------------------------------------------------------------------
        # 1. GENERATE GROUPS AND COMPS FILES (Fixed Logic)
        # ---------------------------------------------------------------------
        echo "Generating Group definitions..." > {log}
        mkdir -p {params.output_dir}/ExpressionInput
        
        group_file="{params.output_dir}/ExpressionInput/groups.{params.task}.txt"
        comp_file="{params.output_dir}/ExpressionInput/comps.{params.task}.txt"
        
        # Clear file first
        > $group_file
        
        for bed_file in {input.junc_beds}; do
            # Extract filename (e.g., CS007_tumor__junction.bed)
            filename=$(basename $bed_file)
            
            # Convert to AltAnalyze expected format (CS007_tumor.bed)
            stream=$(echo $filename | sed 's/__junction.bed/.bed/g')
            
            # EXPLICIT CHECK: Tumor vs Normal
            # We assume the file path or name contains "tumor" or "normal" (case insensitive)
            
            if [[ "$filename" == *"tumor"* ]] || [[ "$filename" == *"Tumor"* ]]; then
                # Group 2 = Tumor (Experiment)
                echo -e "${{stream}}\t2\tTumor" >> $group_file
                echo "Assigned $filename to Group 2 (Tumor)" >> {log}
                
            elif [[ "$filename" == *"normal"* ]] || [[ "$filename" == *"Normal"* ]]; then
                # Group 1 = Normal (Control)
                echo -e "${{stream}}\t1\tNormal" >> $group_file
                echo "Assigned $filename to Group 1 (Normal)" >> {log}
                
            else
                echo "WARNING: Could not determine sample type for $filename. Skipping." >> {log}
            fi
        done

        # Create comps file
        # Compare Group 1 (Normal) vs Group 2 (Tumor)
        echo -e "1\t2" > $comp_file

        # ---------------------------------------------------------------------
        # 2. RUN MULTIPATH-PSI
        # ---------------------------------------------------------------------
        echo "Running AltAnalyze MultiPath-PSI..." >> {log}
        
        # Link database
        if [ -d "$TOOL_DIR/AltDatabase" ]; then rm -rf "$TOOL_DIR/AltDatabase"; fi
        ln -s {params.db_path} "$TOOL_DIR/AltDatabase"
        
        python "$TOOL_DIR/AltAnalyze.py" \
            --species {params.species} \
            --platform {params.platform} \
            --version {params.version} \
            --bedDir {params.bed_dir} \
            --output {params.output_dir} \
            --groupdir "$group_file" \
            --compdir "$comp_file" \
            --expname {params.task} \
            --runGOElite no \
            >> {log} 2>&1

        # ---------------------------------------------------------------------
        # 3. PRUNING
        # ---------------------------------------------------------------------
        echo "Pruning Results..." >> {log}
        
        INPUT_COUNTS="{params.output_dir}/ExpressionInput/counts.{params.task}.txt"
        ANNOTATION="{params.output_dir}/AltResults/AlternativeOutput/{params.species}_{params.platform}_top_alt_junctions-PSI_EventAnnotation.txt"
        OUTPUT_FULL="{params.output_dir}/ExpressionInput/counts.{params.task}.full.txt"
        OUTPUT_PRUNED="{output.counts_pruned_file}"
        OUTPUT_SPECIFIC="{output.counts_pruned_spec_file}"
        
        # Call python script with new arguments
        # Adjust --min_tumor_count and --max_normal_count as needed
        python {input.prune_script} \
            --input "$INPUT_COUNTS" \
            --annotation "$ANNOTATION" \
            --output_full "$OUTPUT_FULL" \
            --output_pruned "$OUTPUT_PRUNED" \
            --output_specific "$OUTPUT_SPECIFIC" \
            --min_tumor_count 10 \
            --max_normal_count 0 \
            >> {log} 2>&1
        """

###########################################################
###########################################################
# ---------------------------------------------------------------------
# 定义所有的剪接类型列表，避免重复书写
# ---------------------------------------------------------------------
SPLICE_TYPES = ["intron-retention", "cassette-exon", "alt-3", "alt-5", "alt-C-term", "altPromoter", "trans-splicing"]

rule Splice_Event_Filter:
    input:
        counts_pruned_spec_file = "rna/altanalyze/output/ExpressionInput/counts.original.specific.txt",
        outFilterbam = "rna/align/{sample}_tumor/Aligned.sortedByCoord.out.bam",
        altsplice_filter_script = workflow.source_path("../scripts/Splice_Event_Filter.py")
    output:
        # 修改点 1: 在 output 中使用 expand 时，如果想保留 {sample} 作为通配符，
        # 需要使用双大括号 {{sample}}。这样 expand 只会展开 splice_type，
        # 而将 {{sample}} 转换为 {sample} 留给 Snakemake 去匹配。
        expand("rna/altanalyze/output/AltAnalyzeFilter/{{sample}}_{splice_type}_Raw.txt",
               splice_type=SPLICE_TYPES)
    params:
        intronJunction = "rna/altanalyze/bed/{sample}_tumor__intronJunction.bed",
        junctionbed =  "rna/altanalyze/bed/{sample}_tumor__junction.bed",
        EventAnnotation = "rna/altanalyze/output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt",
        Gtf_Bed_File_Dir = config["reference"]["Gtf_Bed_File_Dir"],
        outdir = "rna/altanalyze/output/AltAnalyzeFilter"
    resources:
        mem_mb=20000
    log: 
        "logs/altanalyze/AltAnalyze_main.{sample}.log"
    conda:
         "../../envs/bedtools.yaml"
    shell:
        """
        python {input.altsplice_filter_script} \
                        -p {input.counts_pruned_spec_file} \
                        -b {input.outFilterbam} \
                        -i {params.intronJunction} \
                        -j {params.junctionbed} \
                        -e {params.EventAnnotation} \
                        -o {params.outdir} \
                        -s {wildcards.sample} \
                        -g {params.Gtf_Bed_File_Dir} \
                            2>{log}
        """

###########################################################
###########################################################

# 修改点 2: 编写一个输入函数 (Input Function)
# 在这个函数内部，wildcards 对象是可用的。
# 这样我们可以根据当前的 {sample} 动态生成所需的 7 个文件列表。
def get_splice_event_inputs(wildcards):
    return expand("rna/altanalyze/output/AltAnalyzeFilter/{sample}_{splice_type}_Raw.txt",
                  sample=wildcards.sample,
                  splice_type=SPLICE_TYPES)

rule Splice_Event_Specific:
    input:
        # 修改点 3: 将 input 指向上面定义的函数，而不是直接调用 expand
        Splice_Event_Inputs = get_splice_event_inputs,
        counts_pruned_file = "rna/altanalyze/output/ExpressionInput/counts.original.pruned.txt",
        script_path = workflow.source_path("../scripts/Splice_Event_Specific.py")
    output:
        Splice_Event_Specific_Out = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_Specific_High.xlsx"
    params:
        outdir = "rna/altanalyze/output/AltAnalyzeFilter",
        Splice_Event_List_File = "rna/altanalyze/output/AltAnalyzeFilter/Splice_Event_File_Info.txt", # 建议加上 sample 防止并行冲突
        Sample_Info_File = "rna/altanalyze/output/AltAnalyzeFilter/Sample_Info.txt",       # 建议加上 sample 防止并行冲突
        snaf_control = config["reference"]["SNAF_Control"],
        cssw_control = config["reference"]["SRA_Control"],
    resources:
        mem_mb = 20000
    log:
        "logs/altanalyze/alternative_splice_specific.{sample}.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        # 1. 创建输入文件列表
        # 注意: input.Splice_Event_Inputs 是一个列表，直接打印会被空格分开
        # 使用 printf "%s\\n" 可以确保每个文件占一行
        printf "%s\\n" {input.Splice_Event_Inputs} > {params.Splice_Event_List_File} 2> {log}

        # 2. 创建样本信息文件
        echo -e "SampleID\\tBioProject\\tdata_type" > {params.Sample_Info_File}
        echo -e "{wildcards.sample}\\tTargetProject\\trna" >> {params.Sample_Info_File}

        # 3. 运行脚本
        python {input.script_path} \
            -s {params.Splice_Event_List_File} \
            -p {input.counts_pruned_file} \
            -o {params.outdir} \
            -d {params.snaf_control} \
            -c {params.cssw_control} \
            -u {params.Sample_Info_File} \
            -n {wildcards.sample} \
            -a 3 -b 3 \
            >> {log} 2>&1

        # 4. 简单验证
        ls -la {output.Splice_Event_Specific_Out} >> {log} 2>&1
        """

# =============================================================================
# 自动安装外部软件 (TRANSAID & GMST)
# =============================================================================
# =============================================================================
# 0. 自动安装外部软件 (TRANSAID & GMST) - wget 增强版
# =============================================================================
rule install_external_tools:
    output:
        gmst_ok = "resources/software/gmst/gmst.pl",
        transaid_ok = "resources/software/TRANSAID/predict.py",
        transaid_model = "resources/software/TRANSAID/model/TRANSAID_Embedding_batch4_best_model.pth"
    params:
        soft_dir = "resources/software"
    conda:
        "../../envs/translation.yaml"
    shell:
        """
        mkdir -p {params.soft_dir}
        
        # --- 辅助函数 ---
        function git_clone_retry {{
            local repo=$1
            local dir=$2
            if [ -d "$dir" ]; then return 0; fi
            echo "Cloning $repo..."
            # 关闭 SSL 验证以提高国内连接成功率
            git config --global http.sslVerify false
            until git clone --depth 1 "$repo" "$dir"; do
                echo "Clone failed. Retrying in 5 seconds..."
                rm -rf "$dir"
                sleep 5
            done
        }}

        # 1. 安装 GeneMarkS-T (GMST)
        if [ ! -f "{output.gmst_ok}" ]; then
            rm -rf {params.soft_dir}/SQANTI2_tmp
            git_clone_retry https://github.com/Magdoll/SQANTI2.git {params.soft_dir}/SQANTI2_tmp
            
            rm -rf {params.soft_dir}/gmst
            mv {params.soft_dir}/SQANTI2_tmp/utilities/gmst {params.soft_dir}/
            rm -rf {params.soft_dir}/SQANTI2_tmp
            
            chmod +x {params.soft_dir}/gmst/*.pl
            chmod +x {params.soft_dir}/gmst/gmhmmp
            chmod +x {params.soft_dir}/gmst/probuild
        fi

        # 2. 安装 TRANSAID 代码 (仅代码)
        if [ ! -d "{params.soft_dir}/TRANSAID" ]; then
            git_clone_retry https://github.com/wuzengding/TRANSAID.git {params.soft_dir}/TRANSAID
        fi
        
        # 3. 下载 TRANSAID 模型 (使用 wget 直接下载，不依赖 git lfs)
        # 模型文件较大，如果 git clone 下来的是指针文件(1KB左右)，或者文件不存在，则重新下载
        MODEL_FILE="{output.transaid_model}"
        
        if [ ! -f "$MODEL_FILE" ] || [ $(stat -c%s "$MODEL_FILE") -lt 1000000 ]; then
            echo "Downloading TRANSAID model via wget..."
            # 确保目录存在
            mkdir -p $(dirname "$MODEL_FILE")
            
            # 直接从 GitHub Raw/LFS 地址下载
            # 注意：GitHub 的大文件通常可以通过 raw/main/... 获取重定向链接
            wget -O "$MODEL_FILE" "https://github.com/wuzengding/TRANSAID/raw/main/model/TRANSAID_Embedding_batch4_best_model.pth" || \
            wget -O "$MODEL_FILE" "https://media.githubusercontent.com/media/wuzengding/TRANSAID/main/model/TRANSAID_Embedding_batch4_best_model.pth"
            
            # 校验下载是否成功
            if [ ! -f "$MODEL_FILE" ] || [ $(stat -c%s "$MODEL_FILE") -lt 1000000 ]; then
                 echo "ERROR: Model download failed. Please download manually and place at: $MODEL_FILE"
                 exit 1
            fi
        fi
        
        touch {output.transaid_ok}
        """

# =============================================================================
# 1. 翻译预测 (Splice_Event_Trans)
# =============================================================================
rule Splice_Event_Trans:
    input:
        # 上一步的 Excel 结果
        Splice_Event_Specific_Out = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_Specific_High.xlsx",
        script_path = workflow.source_path("../scripts/Splice_Event_Trans.py"),
        # 依赖软件安装
        gmst_tool = "resources/software/gmst/gmst.pl",
        transaid_tool = "resources/software/TRANSAID/predict.py",
        transaid_model = "resources/software/TRANSAID/model/TRANSAID_Embedding_batch4_best_model.pth"
    output:
        Splice_Event_Trans_out = "rna/altanalyze/output/translation/{sample}_Alt_Splice_full_length_trans.txt",
        paper_bed_intersect = "rna/altanalyze/output/translation/{sample}_Alt_Splice_full_length_trans_Intersect_Paper.bed",
        # 标记文件
        trans_done = "rna/altanalyze/output/translation/{sample}_trans.done"
    params:
        outdir = "rna/altanalyze/output/translation",
        sample_name = "{sample}",
        
        # 从 config.yaml 获取的引用数据路径
        Gtf_Bed_File_Dir = config["reference"]["Gtf_Bed_File_Dir"],
        Human_Ref = config["reference"]["genome"],  # 对应 config 中的 genome 路径
        Paper_Bed = config["reference"]["Paper_Bed"], # 需要你在 config 中添加此项
        
        # 软件调用路径 (使用 install rule 下载的相对路径)
        perl = "perl", # 使用 conda 环境中的 perl
        python_transaid = "python", # 使用 conda 环境中的 python
        transaid = "resources/software/TRANSAID/predict.py",
        gmst = "resources/software/gmst/gmst.pl",
        TransDecoder_LongOrfs = "TransDecoder.LongOrfs", # Conda 安装
        TransDecoder_Predict = "TransDecoder.Predict",   # Conda 安装
        model_path = "resources/software/TRANSAID/model/TRANSAID_Embedding_batch4_best_model.pth"
    resources:
        mem_mb = 30000
    threads: 10
    log:
        "logs/altanalyze/splice_trans.{sample}.log"
    conda:
        "../../envs/translation.yaml"
    shell:
        """
        python {input.script_path} \
            -s {input.Splice_Event_Specific_Out} \
            -o {params.outdir} \
            -n {params.sample_name} \
            -g {params.Gtf_Bed_File_Dir} \
            -r {params.Human_Ref} \
            -e {params.perl} \
            -p {params.python_transaid} \
            -t {params.transaid} \
            -m {params.gmst} \
            -l {params.TransDecoder_LongOrfs} \
            -d {params.TransDecoder_Predict} \
            -md {params.model_path} \
            -b {params.Paper_Bed} \
            >> {log} 2>&1
        
        touch {output.trans_done}
        """

# =============================================================================
# 2. 肽段生成 (Splice_Event_Peptide)
# =============================================================================
rule Splice_Event_Peptide:
    input:
        Splice_Event_Trans_out = "rna/altanalyze/output/translation/{sample}_Alt_Splice_full_length_trans.txt",
        paper_bed_intersect = "rna/altanalyze/output/translation/{sample}_Alt_Splice_full_length_trans_Intersect_Paper.bed",
        trans_done = "rna/altanalyze/output/translation/{sample}_trans.done",
        script_path = workflow.source_path("../scripts/Splice_Event_Peptide.py")
    output:
        Splice_Event_Peptide_out = "rna/altanalyze/output/translation/{sample}_Full_Length_Protein_Peptide_info.txt",
        Splice_Event_Peptide_seq_out = "rna/altanalyze/output/translation/{sample}_Full_Length_Protein_Peptide_seq.faa"
    params:
        outdir = "rna/altanalyze/output/translation",
        sample_name = "{sample}",
        Gtf_Bed_File_Dir = config["reference"]["Gtf_Bed_File_Dir"],
        Paper_Bed = config["reference"]["Paper_Bed"]
    resources:
        mem_mb = 10000
    log:
        "logs/altanalyze/splice_peptide.{sample}.log"
    conda:
        "../../envs/translation.yaml"
    shell:
        """
        python {input.script_path} \
            -s {input.Splice_Event_Trans_out} \
            -o {params.outdir} \
            -n {params.sample_name} \
            -g {params.Gtf_Bed_File_Dir} \
            -b {input.paper_bed_intersect} \
            >> {log} 2>&1
        """
# =============================================================================
# Rule: Splice_Event_Final_Report
# 功能: 整合特异性分析结果与多肽预测结果，生成最终报告
# =============================================================================
rule Splice_Event_Final_Report:
    input:
        # 输入1: 特异性筛选 Excel (来自 Splice_Event_Specific)
        specific_xlsx = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_Specific_High.xlsx",
        # 输入2: 肽段信息 (来自 Splice_Event_Peptide)
        peptide_info = "rna/altanalyze/output/translation/{sample}_Full_Length_Protein_Peptide_info.txt",
        script_path = workflow.source_path("../scripts/Splice_Event_Final_Report.py")
    output:
        final_xlsx = "rna/altanalyze/output/translation/{sample}_Final_Neoantigen_Report.xlsx"
    params:
        # 无需额外参数，直接传递文件路径
    resources:
        mem_mb = 4000
    log:
        "logs/altanalyze/final_report.{sample}.log"
    conda:
        "../../envs/translation.yaml" # 只需要 pandas 和 openpyxl
    shell:
        """
        python {input.script_path} \
            -s {input.specific_xlsx} \
            -p {input.peptide_info} \
            -o {output.final_xlsx} \
            >> {log} 2>&1
        """