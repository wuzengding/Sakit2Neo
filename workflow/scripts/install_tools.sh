#!/bin/bash

#确保以root身份运行
if [ "$(id -u)" -ne 0 ]; then
   echo "此脚本必须以root身份运行" >&2
   echo "请尝试使用: sudo $0" >&2
   exit 1
fi

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"
#!/bin/bash

# 创建日志文件
LOG_FILE="tools_installation.log"
echo "开始安装生物信息学工具 $(date)" > $LOG_FILE

# 创建安装目录
mkdir -p /opt/tools
mkdir -p /opt/tools/bin

# 设置PATH
echo 'export PATH=/opt/tools/bin:$PATH' >> ~/.bashrc

# 函数：记录日志消息
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $LOG_FILE
}

# 安装基础软件包
log "安装基础软件包..."
apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    python3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

if [ $? -eq 0 ]; then
    log "基础软件包安装成功"
else
    log "基础软件包安装失败"
    exit 1
fi

# 设置Python别名
ln -sf /usr/bin/python3 /usr/bin/python
ln -sf /usr/bin/pip3 /usr/bin/pip

# 安装Python包
log "安装Python必要包..."
pip install biopython pandas matplotlib numpy scipy seaborn pysam

# 安装Samtools
log "安装Samtools..."
cd /opt/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/opt/tools
make
make install
ln -sf /opt/tools/bin/samtools /opt/tools/bin/
if [ -f /opt/tools/bin/samtools ]; then
    log "Samtools安装成功: $(samtools --version | head -n 1)"
else
    log "Samtools安装失败"
fi

# 安装BWA
log "安装BWA..."
cd /opt/tools
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa /opt/tools/bin/
if [ -f /opt/tools/bin/bwa ]; then
    log "BWA安装成功: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
else
    log "BWA安装失败"
fi

# 安装GATK
log "安装GATK..."
cd /opt/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
ln -sf /opt/tools/gatk-4.2.6.1/gatk /opt/tools/bin/
if [ -f /opt/tools/bin/gatk ]; then
    log "GATK安装成功: $(gatk --version | head -n 1)"
else
    log "GATK安装失败"
fi

# 安装VarScan
log "安装VarScan..."
cd /opt/tools
wget https://github.com/dkoboldt/varscan/releases/download/2.4.4/VarScan.v2.4.4.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/VarScan.v2.4.4.jar "$@"' > /opt/tools/bin/varscan
chmod +x /opt/tools/bin/varscan
if [ -f /opt/tools/bin/varscan ]; then
    log "VarScan安装成功: $(varscan 2>&1 | grep VarScan | head -n 1)"
else
    log "VarScan安装失败"
fi

# 安装Strelka2
log "安装Strelka2..."
cd /opt/tools
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
ln -sf /opt/tools/strelka-2.9.10.centos6_x86_64/bin/* /opt/tools/bin/
if [ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ]; then
    log "Strelka2安装成功"
else
    log "Strelka2安装失败"
fi

# 安装bcftools
log "安装bcftools..."
cd /opt/tools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjf bcftools-1.16.tar.bz2
cd bcftools-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/bcftools ]; then
    log "bcftools安装成功: $(bcftools --version | head -n 1)"
else
    log "bcftools安装失败"
fi

# 安装htslib
log "安装htslib..."
cd /opt/tools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/opt/tools
make
make install
if [ -f /opt/tools/bin/tabix ]; then
    log "htslib安装成功: $(tabix --version | head -n 1)"
else
    log "htslib安装失败"
fi

# 安装FastQC
log "安装FastQC..."
cd /opt/tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 /opt/tools/FastQC/fastqc
ln -sf /opt/tools/FastQC/fastqc /opt/tools/bin/
if [ -f /opt/tools/bin/fastqc ]; then
    log "FastQC安装成功: $(fastqc -v)"
else
    log "FastQC安装失败"
fi

# 安装Trimmomatic
log "安装Trimmomatic..."
cd /opt/tools
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/Trimmomatic-0.39/trimmomatic-0.39.jar "$@"' > /opt/tools/bin/trimmomatic
chmod +x /opt/tools/bin/trimmomatic
if [ -f /opt/tools/bin/trimmomatic ]; then
    log "Trimmomatic安装成功"
else
    log "Trimmomatic安装失败"
fi

# 安装Picard
log "安装Picard..."
cd /opt/tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar
echo -e '#!/bin/bash\njava -jar /opt/tools/picard.jar "$@"' > /opt/tools/bin/picard
chmod +x /opt/tools/bin/picard
if [ -f /opt/tools/bin/picard ]; then
    log "Picard安装成功"
else
    log "Picard安装失败"
fi

# 安装SnpEff
log "安装SnpEff..."
cd /opt/tools
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/snpEff.jar "$@"' > /opt/tools/bin/snpeff
chmod +x /opt/tools/bin/snpeff
echo -e '#!/bin/bash\njava -jar /opt/tools/snpEff/SnpSift.jar "$@"' > /opt/tools/bin/snpsift
chmod +x /opt/tools/bin/snpsift
if [ -f /opt/tools/bin/snpeff ]; then
    log "SnpEff安装成功"
else
    log "SnpEff安装失败"
fi

# 清理下载文件
log "清理下载文件..."
cd /opt/tools
rm -f *.tar.bz2 *.zip

# 设置权限
chmod -R 755 /opt/tools/bin

# 创建一个测试脚本
cat > /opt/tools/test_tools.sh << 'EOF'
#!/bin/bash

echo "测试生物信息学工具安装:"
echo "---------------------"
echo "Samtools: $(samtools --version | head -n 1)"
echo "BWA: $(bwa 2>&1 | grep Version | cut -f 2 -d ' ')"
echo "GATK: $(gatk --version | head -n 1)"
echo "BCFTools: $(bcftools --version | head -n 1)"
echo "FastQC: $(fastqc -v)"
echo "Python: $(python --version)"
echo "VarScan存在: $([ -f /opt/tools/bin/varscan ] && echo '是' || echo '否')"
echo "Strelka2存在: $([ -f /opt/tools/bin/configureStrelkaSomaticWorkflow.py ] && echo '是' || echo '否')"
echo "Trimmomatic存在: $([ -f /opt/tools/bin/trimmomatic ] && echo '是' || echo '否')"
echo "Picard存在: $([ -f /opt/tools/bin/picard ] && echo '是' || echo '否')"
echo "SnpEff存在: $([ -f /opt/tools/bin/snpeff ] && echo '是' || echo '否')"
echo "---------------------"
EOF

chmod +x /opt/tools/test_tools.sh

# 完成安装
log "所有工具安装完成!"
log "请运行以下命令重新加载环境变量: source ~/.bashrc"
log "然后运行测试脚本检查安装: /opt/tools/test_tools.sh"

echo "安装完成，请执行 'source ~/.bashrc' 更新环境变量，然后执行 '/opt/tools/test_tools.sh' 测试工具安装"

