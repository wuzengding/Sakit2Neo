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