#在docker SakitNeo容器中执行以下代码#
snakemake -s Snakefile --cores 40 --configfile ../config/config.yaml --use-conda --conda-prefix /home/work/SakitNeo_envs/ --rerun-incomplete
