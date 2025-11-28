#在docker SakitNeo容器中执行以下代码#
snakemake -s /home/bioinfo/05.pipeline_dev/SakitNeo/workflow/Snakefile --cores 40 --configfile /home/bioinfo/09.data_CCS_sdfyy/05.Result.Sakit2Neo/CS007/config_CS007.yaml --use-conda --conda-prefix /home/work/SakitNeo_envs/ --rerun-incomplete --keep-going
