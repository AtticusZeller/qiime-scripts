#!/bin/bash
# ref: https://mp.weixin.qq.com/s/KgNzy_TWcQjO4firXeadZQ?scene=1&click_id=1
# =======================================================
# 脚本名称: 01_setup_resources.sh
# 功能: 安装 QIIME2 环境, Aspera, 并构建 SILVA 138.1 分类器
# =======================================================

set -e # 遇到错误立即停止

# --- 1. 配置参数 ---
DB_DIR="$PWD/db"
THREADS=16

mkdir -p $DB_DIR

echo "[INFO] === 第一步：检查并安装 Conda 环境 ==="
# 检查是否已存在名为 qiime2-amplicon-2024.5 的环境
if conda info --envs | grep -q "qiime2-amplicon-2024.5"; then
    echo "[INFO] 环境 qiime2-amplicon-2024.5 已存在，跳过安装。"
else
    echo "[INFO] 正在创建 QIIME 2 环境..."
    conda update -n base -c defaults conda -y
    conda env create -n qiime2-amplicon-2024.5 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml
fi

# 激活环境 (在Shell脚本中激活conda需要特殊处理)
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.5

echo "[INFO] === 第二步：检查并安装 Aspera Connect ==="
export PATH=~/.aspera/connect/bin:$PATH
if ! command -v ascp &> /dev/null; then
    echo "[INFO] 未检测到 ascp，正在安装..."
    wget -qO- https://download.asperasoft.com/download/sw/connect/3.9.8/ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz | tar xvz
    ./ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.sh
    # 将路径永久写入配置 (可选)
    # echo 'export PATH=~/.aspera/connect/bin:$PATH' >> ~/.bashrc
else
    echo "[INFO] Aspera 已安装。"
fi

echo "[INFO] === 第三步：构建 SILVA 138.1 分类器 (RESCRIPt 流程) ==="
cd $DB_DIR

# 3.1 下载原始文件 (如果文件不存在则下载)
echo "[LOG] 正在检查/下载 SILVA 原始文件...关掉代理以防下载失败"
[ ! -f "tax_slv_ssu_138.1.txt" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz && gunzip tax_slv_ssu_138.1.txt.gz
[ ! -f "taxmap_slv_ssu_ref_nr_138.1.txt" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz && gunzip taxmap_slv_ssu_ref_nr_138.1.txt.gz
[ ! -f "tax_slv_ssu_138.1.tre" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.tre.gz && gunzip tax_slv_ssu_138.1.tre.gz
[ ! -f "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz && gunzip SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz

# 3.2 导入 QIIME2
echo "[LOG] 导入数据到 QIIME2..."
qiime tools import --type 'FeatureData[SILVATaxonomy]' --input-path tax_slv_ssu_138.1.txt --output-path taxranks-silva-138.1-ssu-nr99.qza
qiime tools import --type 'FeatureData[SILVATaxidMap]' --input-path taxmap_slv_ssu_ref_nr_138.1.txt --output-path taxmap-silva-138.1-ssu-nr99.qza
qiime tools import --type 'Phylogeny[Rooted]' --input-path tax_slv_ssu_138.1.tre --output-path taxtree-silva-138.1-nr99.qza
qiime tools import --type 'FeatureData[RNASequence]' --input-path SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta --output-path silva-138.1-ssu-nr99-rna-seqs.qza

# 3.3 转换与清洗
echo "[LOG] RESCRIPt 转换与清洗..."
qiime rescript reverse-transcribe --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza
qiime rescript parse-silva-taxonomy --i-taxonomy-tree taxtree-silva-138.1-nr99.qza --i-taxonomy-map taxmap-silva-138.1-ssu-nr99.qza --i-taxonomy-ranks taxranks-silva-138.1-ssu-nr99.qza --o-taxonomy silva-138.1-ssu-nr99-tax.qza

qiime rescript cull-seqs --i-sequences silva-138.1-ssu-nr99-seqs.qza --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.1-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza

qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.1-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza

# 3.4 提取特定引物区段 (338F-806R)
echo "[LOG] 提取 338F-806R 片段..."
qiime feature-classifier extract-reads \
    --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer ACTCCTACGGGAGGCAGCAG \
    --p-r-primer GGACTACHVGGGTWTCTAAT \
    --p-n-jobs $THREADS \
    --p-read-orientation 'forward' \
    --o-reads silva-138.1-ssu-nr99-seqs-338f-806r.qza

qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-338f-806r.qza \
    --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-338f-806r-uniq.qza \
    --o-dereplicated-taxa  silva-138.1-ssu-nr99-tax-338f-806r-derep-uniq.qza \
    --p-threads $THREADS

# 3.5 训练分类器
echo "[LOG] 训练分类器 (这可能需要一些时间)..."
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.1-ssu-nr99-seqs-338f-806r-uniq.qza \
    --i-reference-taxonomy silva-138.1-ssu-nr99-tax-338f-806r-derep-uniq.qza \
    --o-classifier silva-138.1-ssu-nr99-338f-806r-classifier.qza

echo "[SUCCESS] 资源准备完成！分类器位置: $DB_DIR/silva-138.1-ssu-nr99-338f-806r-classifier.qza"