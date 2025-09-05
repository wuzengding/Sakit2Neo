import os
import yaml
import datetime
import subprocess
from snakemake.utils import min_version

# 确保最低Snakemake版本
min_version("6.0.0")

# 处理全局变量
try:
    REFERENCES_TO_INSTALL
    SOFTWARE_TO_INSTALL
except NameError:
    # 如果变量未定义，重新计算
    # 加载环境配置
    env_config_path = "../../config/env_config.yaml"
    with open(env_config_path, 'r') as f:
        env_config = yaml.safe_load(f)
        
    # 导入检测函数
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))
    from scripts.check_env import check_references_and_software
    
    # 提取路径信息
    resources_base_dir = env_config.get("resources_base_dir", "/home/bioinfo/05.pipeline_dev/SakitNeo/resources")
    ref_paths = env_config.get("reference_paths", {})
    ref_base = ref_paths.get("base_dir", os.path.join(resources_base_dir, "references"))
    genomes_dir = ref_paths.get("genomes_dir", os.path.join(ref_base, "genomes"))
    variants_dir = ref_paths.get("variants_dir", os.path.join(ref_base, "variants"))
    annotation_dir = ref_paths.get("annotation_dir", os.path.join(ref_base, "annotation"))
    targets_dir = ref_paths.get("targets_dir", os.path.join(ref_base, "targets"))
    
    # 重新计算需要安装的组件
    REFERENCES_TO_INSTALL, SOFTWARE_TO_INSTALL = check_references_and_software(
        env_config,
        genomes_dir=genomes_dir,
        variants_dir=variants_dir,
        annotation_dir=annotation_dir,
        targets_dir=targets_dir
    )

# 加载环境配置
env_config_path = "../../config/env_config.yaml"
with open(env_config_path, 'r') as f:
    env_config = yaml.safe_load(f)

# 提取配置信息
RESOURCES_BASE_DIR = env_config.get("resources_base_dir", "/home/bioinfo/05.pipeline_dev/SakitNeo/resources")
USE_CONDA = env_config.get("use_conda", True)

# 各种路径
REF_PATHS = env_config.get("reference_paths", {})
REF_BASE_DIR = REF_PATHS.get("base_dir", os.path.join(RESOURCES_BASE_DIR, "references"))
GENOMES_DIR = REF_PATHS.get("genomes_dir", os.path.join(REF_BASE_DIR, "genomes"))
ANNOTATION_DIR = REF_PATHS.get("annotation_dir", os.path.join(REF_BASE_DIR, "annotation"))
VARIANTS_DIR = REF_PATHS.get("variants_dir", os.path.join(REF_BASE_DIR, "variants"))
TARGETS_DIR = REF_PATHS.get("targets_dir", os.path.join(REF_BASE_DIR, "targets"))

SW_PATHS = env_config.get("software_paths", {})
CONDA_BASE_DIR = SW_PATHS.get("conda_base_dir", os.path.join(RESOURCES_BASE_DIR, "conda"))
TOOLS_DIR = SW_PATHS.get("tools_dir", os.path.join(RESOURCES_BASE_DIR, "tools"))
BIN_DIR = SW_PATHS.get("bin_dir", os.path.join(TOOLS_DIR, "bin"))

# 日志目录
LOG_DIR = os.path.join(RESOURCES_BASE_DIR, "logs")

# 创建所需目录
def create_required_dirs():
    """创建所有必需的目录"""
    for directory in [RESOURCES_BASE_DIR, REF_BASE_DIR, GENOMES_DIR, ANNOTATION_DIR, 
                      VARIANTS_DIR, TARGETS_DIR, TOOLS_DIR, BIN_DIR, LOG_DIR]:
        os.makedirs(directory, exist_ok=True)
    
    if USE_CONDA:
        os.makedirs(CONDA_BASE_DIR, exist_ok=True)

# 环境配置标记文件
ENV_SETUP_DONE = os.path.join(RESOURCES_BASE_DIR, "env_setup_complete.txt")

# 创建一个标记文件，表示环境已经配置完成
rule env_setup:
    output:
        ENV_SETUP_DONE
    log:
        os.path.join(LOG_DIR, "env_setup.log")
    threads: 4
    resources:
        mem_mb=4000,
        time=120
    run:
        # 创建必需的目录
        create_required_dirs()
        
        # 安装conda环境
        if "conda_env" in SOFTWARE_TO_INSTALL:
            conda_env_name = env_config.get("software", {}).get("conda_env_name", "bioinfo_env")
            channels = " ".join([f"-c {c}" for c in env_config.get("software", {}).get("conda_channels", [])])
            
            packages = []
            for pkg_name, pkg_config in env_config.get("software", {}).get("packages", {}).items():
                if pkg_config.get("conda", True) and pkg_config.get("enabled", True):
                    version = pkg_config.get("version", "")
                    if version:
                        packages.append(f"{pkg_name}={version}")
                    else:
                        packages.append(pkg_name)
            
            pkg_str = " ".join(packages)
            
            shell(f"""
                # 创建conda环境
                export CONDA_PKGS_DIRS={CONDA_BASE_DIR}/pkgs
                export CONDA_ENVS_PATH={CONDA_BASE_DIR}/envs
                conda create -y -p {CONDA_BASE_DIR}/envs/{conda_env_name} {channels} {pkg_str} &> {log}
            """)
        
        # 安装非conda软件
        for sw_name, sw_config in env_config.get("software", {}).get("packages", {}).items():
            if not sw_config.get("conda", True) and sw_config.get("enabled", True) and sw_name in SOFTWARE_TO_INSTALL:
                install_cmd = sw_config.get("install_cmd", "")
                if install_cmd:
                    install_cmd = install_cmd.format(
                        software_paths=SW_PATHS,
                        reference_paths=REF_PATHS
                    )
                    shell(f"""
                        # 安装{sw_name}
                        {install_cmd} &>> {log}
                    """)
        
        # 下载参考基因组
        for ref_name, ref_config in env_config.get("references", {}).items():
            if not ref_config.get("enabled", True) or ref_name not in REFERENCES_TO_INSTALL:
                continue
                
            source = ref_config.get("source", "")
            dest = ref_config.get("destination", "")
            
            # 替换路径中的变量
            dest = dest.format(
                genomes_dir=GENOMES_DIR,
                variants_dir=VARIANTS_DIR,
                annotation_dir=ANNOTATION_DIR,
                targets_dir=TARGETS_DIR
            )
            
            # 创建目标所在目录
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            
            # 处理需要认证的资源
            if ref_config.get("requires_authentication", False):
                shell(f"""
                    echo "警告: {ref_name} 需要认证，请手动下载到 {dest}"
                    echo "下载URL: {source}"
                    echo "认证URL: {ref_config.get('auth_url', '')}"
                """)
                continue
            
            # 下载文件
            shell(f"wget -O {dest} {source} &>> {log}")
            
            # 解压缩
            if ref_config.get("decompress", False) and dest.endswith('.gz'):
                shell(f"gunzip -f {dest} &>> {log}")
            
            # 处理命令
            process_cmd = ref_config.get("process_cmd", "")
            if process_cmd:
                process_cmd = process_cmd.format(
                    genomes_dir=GENOMES_DIR,
                    variants_dir=VARIANTS_DIR, 
                    annotation_dir=ANNOTATION_DIR,
                    targets_dir=TARGETS_DIR
                )
                shell(f"{process_cmd} &>> {log}")
            
            # 创建索引
            for idx_type in ref_config.get("indexes", []):
                if idx_type == "bwa":
                    shell(f"bwa index {dest} &>> {log}")
                elif idx_type == "samtools":
                    shell(f"samtools faidx {dest} &>> {log}")
                elif idx_type == "picard":
                    shell(f"picard CreateSequenceDictionary R={dest} O={os.path.splitext(dest)[0]}.dict &>> {log}")
                elif idx_type == "tabix":
                    shell(f"tabix -p vcf {dest} &>> {log}")
        
        # 执行后处理命令
        post_cmds = env_config.get("post_setup_commands", "")
        if post_cmds:
            post_cmds = post_cmds.format(
                resources_base_dir=RESOURCES_BASE_DIR,
                software_paths=SW_PATHS,
                reference_paths=REF_PATHS
            )
            shell(f"{post_cmds} &>> {log}")
        
        # 创建标记文件
        with open(output[0], "w") as f:
            f.write(f"环境配置完成于 {datetime.datetime.now()}\n")
            f.write(f"安装的参考基因组: {', '.join(REFERENCES_TO_INSTALL) if REFERENCES_TO_INSTALL else '无'}\n")
            f.write(f"安装的软件: {', '.join(SOFTWARE_TO_INSTALL) if SOFTWARE_TO_INSTALL else '无'}\n")
