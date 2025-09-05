# scripts/check_env.py
import os
import subprocess
import yaml

def check_references_and_software(env_config, genomes_dir, variants_dir, annotation_dir, targets_dir):
    """
    检查参考基因组文件和软件是否需要安装
    
    参数:
    env_config: 环境配置字典
    genomes_dir: 基因组目录路径
    variants_dir: 变异数据库目录路径
    annotation_dir: 注释数据目录路径
    targets_dir: 靶向数据目录路径
    
    返回:
    需要安装的参考基因组列表和需要安装的软件列表
    """
    references_to_install = []
    software_to_install = []
    
    # 检查参考基因组
    for ref_name, ref_config in env_config.get("references", {}).items():
        if not ref_config.get("enabled", True):
            continue
            
        dest = ref_config.get("destination", "")
        # 替换路径中的变量
        dest = dest.format(
            genomes_dir=genomes_dir,
            variants_dir=variants_dir,
            annotation_dir=annotation_dir,
            targets_dir=targets_dir
        )
        
        if not os.path.exists(dest):
            references_to_install.append(ref_name)
    
    # 检查软件
    use_conda = env_config.get("use_conda", True)
    if use_conda:
        # 如果使用conda，检查环境是否存在
        conda_env_name = env_config.get("software", {}).get("conda_env_name", "bioinfo_env")
        try:
            result = subprocess.run(['conda', 'env', 'list'], 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, 
                                   text=True)
            if conda_env_name not in result.stdout:
                software_to_install.append("conda_env")
        except:
            software_to_install.append("conda_env")
    else:
        # 如果不使用conda，检查每个工具是否安装
        for sw_name, sw_config in env_config.get("software", {}).get("packages", {}).items():
            if not sw_config.get("enabled", True):
                continue
                
            check_cmd = sw_config.get("check_cmd", "")
            if not check_cmd:
                continue
                
            try:
                result = subprocess.run(check_cmd, shell=True, 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE,
                                       text=True)
                if result.returncode != 0:
                    software_to_install.append(sw_name)
            except:
                software_to_install.append(sw_name)
    
    return references_to_install, software_to_install
