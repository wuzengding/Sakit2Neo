# Helper functions (Keep these as they were)
def get_rna_samples():
    return samples[samples['data_type'] == 'rna']

def get_all_bed_files(file_type="junction"):
    files = []
    for index, row in get_rna_samples().iterrows():
        s = row["sample_id"]
        t = row["sampletype"]
        files.append(f"rna/altanalyze/bed/{s}_{t}__{file_type}.bed")
    return files

# 1. BAM TO BED (Unchanged, runs per sample)
rule bam_to_bed:
    input:
        bam = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam",
        bai = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam.bai"
    output:
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
        TOOL_DIR="$CONDA_PREFIX/share/AltAnalyze"
        if [ ! -f "$TOOL_DIR/AltAnalyze.py" ]; then
            mkdir -p "$CONDA_PREFIX/share"
            git clone https://github.com/nsalomonis/altanalyze.git "$TOOL_DIR" >> {log} 2>&1
        fi
        
        TEMP_BAM="{params.bed_dir}/{wildcards.sample}_{wildcards.sampletype}.bam"
        mkdir -p {params.bed_dir}
        ln -sf $(readlink -f {input.bam}) "$TEMP_BAM"
        ln -sf $(readlink -f {input.bai}) "$TEMP_BAM.bai"
        
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

# 2. ALTANALYZE RUN (Fixed: Runs ONCE for the whole cohort)
rule AltAnalyze_run:
    input:
        junc_beds = get_all_bed_files("junction"),
        exon_beds = get_all_bed_files("intronJunction")
    output:
        # CHANGED: Output is now a single combined file, not per-sample files.
        counts_pruned_file = "rna/altanalyze/output/ExpressionInput/counts.original.pruned.txt",
        counts_full_file = "rna/altanalyze/output/ExpressionInput/counts.original.full.txt"
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
        prune_script = workflow.source_path("../scripts/prune_altanalyze.py")
    conda:
        "../../envs/altanalyze.yaml"
    shell:
        """
        TOOL_DIR="$CONDA_PREFIX/share/AltAnalyze"
        
        echo "Generating Group definitions..." > {log}
        mkdir -p {params.output_dir}/ExpressionInput
        
        group_file="{params.output_dir}/ExpressionInput/groups.{params.task}.txt"
        comp_file="{params.output_dir}/ExpressionInput/comps.{params.task}.txt"
        
        > $group_file
        
        # Iterate over inputs to build the groups file
        for bed_file in {input.junc_beds}; do
            filename=$(basename $bed_file)
            stream=$(echo $filename | sed 's/__junction.bed/.bed/g')
            
            if [[ "$filename" == *"tumor"* ]] || [[ "$filename" == *"Tumor"* ]]; then
                echo -e "${{stream}}\t2\tTumor" >> $group_file
            elif [[ "$filename" == *"normal"* ]] || [[ "$filename" == *"Normal"* ]]; then
                echo -e "${{stream}}\t1\tNormal" >> $group_file
            fi
        done

        echo -e "1\t2" > $comp_file

        echo "Running AltAnalyze MultiPath-PSI..." >> {log}
        
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

        echo "Pruning Results..." >> {log}
        
        INPUT_COUNTS="{params.output_dir}/ExpressionInput/counts.{params.task}.txt"
        ANNOTATION="{params.output_dir}/AltResults/AlternativeOutput/{params.species}_{params.platform}_top_alt_junctions-PSI_EventAnnotation.txt"
        
        python {params.prune_script} \
            --input "$INPUT_COUNTS" \
            --annotation "$ANNOTATION" \
            --output_full "{output.counts_full_file}" \
            --output_pruned "{output.counts_pruned_file}" \
            >> {log} 2>&1
        """

# 3. SPLICE EVENT FILTER (Fixed: Accepts combined matrix, filters for specific sample)
rule Splice_Event_Filter:
    input:
        # CHANGED: Input is the combined matrix from the previous rule
        prunedcount = "rna/altanalyze/output/ExpressionInput/counts.original.pruned.txt",
        outFilterbam = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam",
    output:
        Splice_Filter_Out = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_intron-retention_Raw.txt"
    params:
        outAltAnalyze_dir = "rna/altanalyze/output",
        # Note: These bed paths must match what bam_to_bed outputted
        intronJunction = "rna/altanalyze/bed/{sample}_tumor__intronJunction.bed",
        junctionbed =  "rna/altanalyze/bed/{sample}_tumor__junction.bed",
        EventAnnotation = "rna/altanalyze/output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt",
        Gtf_Bed_File_Dir = config["reference"]["Gtf_Bed_File_Dir"],
        outdir = "rna/altanalyze/output/AltAnalyzeFilter",
        altsplice_filter_script = workflow.source_path("../scripts/Splice_Event_Filter.py")
    resources:
        mem_mb=20000
    log: 
        splice_filter_log = "logs/altanalyze/filter_{sample}.log"
    conda:
         "../../envs/altanalyze.yaml"
    shell:
        """
        # Ensure the filter script knows to look for this specific sample in the big matrix
        python {params.altsplice_filter_script} \
            -p {input.prunedcount} \
            -b {input.outFilterbam} \
            -i {params.intronJunction} \
            -j {params.junctionbed} \
            -e {params.EventAnnotation} \
            -o {params.outdir} \
            -s {wildcards.sample} \
            -g {params.Gtf_Bed_File_Dir} \
            2> {log.splice_filter_log}
        """

# 4. SPLICE EVENT SPECIFIC (Fixed: Correct Output/Directory handling)
rule Splice_Event_Specific:
    input:
        Splice_Event_Input = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_intron-retention_Raw.txt",
        # Use the combined matrix here too
        prunedcount = "rna/altanalyze/output/ExpressionInput/counts.original.pruned.txt"
    output:
        # Main result file
        Splice_Event_Specific_Out = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_specific/intron-retention_Specific_High.txt",
        # Using directory to capture all other outputs in that folder
        out_dir_marker = directory("rna/altanalyze/output/AltAnalyzeFilter/{sample}_specific")
    params:
        # Output to a sample-specific folder to prevent overwriting
        outdir = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_specific",
        Splice_Event_List_File = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_specific/Splice_Event_File_Info.txt",
        Sample_Info_File = "rna/altanalyze/output/AltAnalyzeFilter/{sample}_specific/Sample_Info.txt",
        snaf_control = config["reference"]["SNAF_Control"],
        cssw_control = config["reference"]["SRA_Control"],
        script_path = workflow.source_path("../scripts/Splice_Event_Specific.py")
    resources:
        mem_mb = 20000
    log:
        "logs/altanalyze/alternative_splice_specific_{sample}.log"
    conda:
        "../../envs/python.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        # 1. Create file list
        echo "{input.Splice_Event_Input}" > {params.Splice_Event_List_File} 2> {log}

        # 2. Create dummy info file
        echo -e "SampleID\tBioProject\tdata_type" > {params.Sample_Info_File}
        echo -e "{wildcards.sample}\tTargetProject\trna" >> {params.Sample_Info_File}

        # 3. Run script
        python {params.script_path} \
            -s {params.Splice_Event_List_File} \
            -p {input.prunedcount} \
            -o {params.outdir} \
            -d {params.snaf_control} \
            -c {params.cssw_control} \
            -u {params.Sample_Info_File} \
            -a 3 -b 3 \
            >> {log} 2>&1
        """