#!/bin/bash
# ref: https://mp.weixin.qq.com/s/KgNzy_TWcQjO4firXeadZQ?scene=1&click_id=1
# =======================================================
# 脚本名称: 01_setup_resources.sh
# 功能: 安装环境 + 准备 SILVA 138.1 公共资源 + 构建多区域分类器
# =======================================================

set -e 

# --- 1. 配置参数 ---
DB_DIR="$PWD/database"
THREADS=16

mkdir -p $DB_DIR

# --- 2. 环境与 Aspera (保持不变) ---
source ~/miniconda3/etc/profile.d/conda.sh
echo "[INFO] === 第一步：检查并安装 Conda 环境 ==="
if conda info --envs | grep -q "qiime2-amplicon-2024.5"; then
    echo "[INFO] 环境 qiime2-amplicon-2024.5 已存在，跳过安装。"
else
    echo "[INFO] 正在创建 QIIME 2 环境..."
    conda update -n base -c defaults conda -y
    conda env create -n qiime2-amplicon-2024.5 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml
fi

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

# --- 3. 准备 SILVA 公共资源 (只运行一次) ---
echo "[INFO] === 第三步：准备 SILVA 138.1 基础全长数据 (公共部分) ==="
cd $DB_DIR

# 3.1 下载原始文件 (如果文件不存在则下载)
echo "[LOG] 正在检查/下载 SILVA 原始文件...关掉代理以防下载失败"
[ ! -f "tax_slv_ssu_138.1.txt" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz && gunzip tax_slv_ssu_138.1.txt.gz
[ ! -f "taxmap_slv_ssu_ref_nr_138.1.txt" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz && gunzip taxmap_slv_ssu_ref_nr_138.1.txt.gz
[ ! -f "tax_slv_ssu_138.1.tre" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.tre.gz && gunzip tax_slv_ssu_138.1.tre.gz
[ ! -f "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta" ] && wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz && gunzip SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz

# 3.2 导入 & 3.3 清洗 (这部分最耗时，加一个检测跳过)
FULL_DEREP_SEQS="silva-138.1-ssu-nr99-seqs-derep-uniq.qza"
FULL_DEREP_TAXA="silva-138.1-ssu-nr99-tax-derep-uniq.qza"

if [ -f "$FULL_DEREP_SEQS" ] && [ -f "$FULL_DEREP_TAXA" ]; then
    echo "[SKIP] 全长清洗后序列已存在，跳过 3.2 和 3.3 步骤..."
else
    echo "[LOG] 导入数据并执行 RESCRIPt 清洗流程..."
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
        --o-dereplicated-sequences "$FULL_DEREP_SEQS" \
        --o-dereplicated-taxa "$FULL_DEREP_TAXA"
fi


# --- 4. 构建特定区域分类器 (封装函数) ---
# 功能：传入引物名称和序列，自动生成对应的中间文件和分类器
build_classifier() {
    local p_name=$1  # 名字，例如 338f-806r
    local f_seq=$2   # 前向引物序列
    local r_seq=$3   # 反向引物序列

    local output_classifier="silva-138.1-ssu-nr99-${p_name}-classifier.qza"
    
    echo "--------------------------------------------------------"
    echo "[LOG] 开始构建分类器: ${p_name}"
    echo "      Forward: ${f_seq}"
    echo "      Reverse: ${r_seq}"
    
    if [ -f "$output_classifier" ]; then
        echo "[SKIP] 分类器 ${output_classifier} 已存在，跳过。"
        return
    fi

    # 4.1 提取
    echo "[LOG] Extracting reads for ${p_name}..."
    qiime feature-classifier extract-reads \
        --i-sequences "$FULL_DEREP_SEQS" \
        --p-f-primer "$f_seq" \
        --p-r-primer "$r_seq" \
        --p-n-jobs $THREADS \
        --p-read-orientation 'forward' \
        --o-reads "silva-138.1-ssu-nr99-seqs-${p_name}.qza"

    # 4.2 去重 (针对该片段)
    echo "[LOG] Dereplicating for ${p_name}..."
    qiime rescript dereplicate \
        --i-sequences "silva-138.1-ssu-nr99-seqs-${p_name}.qza" \
        --i-taxa "$FULL_DEREP_TAXA" \
        --p-mode 'uniq' \
        --o-dereplicated-sequences "silva-138.1-ssu-nr99-seqs-${p_name}-uniq.qza" \
        --o-dereplicated-taxa  "silva-138.1-ssu-nr99-tax-${p_name}-derep-uniq.qza" \
        --p-threads $THREADS

    # 4.3 训练
    echo "[LOG] Training classifier for ${p_name}..."
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "silva-138.1-ssu-nr99-seqs-${p_name}-uniq.qza" \
        --i-reference-taxonomy "silva-138.1-ssu-nr99-tax-${p_name}-derep-uniq.qza" \
        --o-classifier "$output_classifier"
        
    echo "[SUCCESS] ${p_name} 分类器构建完成！"
}


# === 5. 执行具体构建任务 ===

# V3-V4 (338F - 806R)
# 之前的: ACTCCTACGGGAGGCAGCAG / GGACTACHVGGGTWTCTAAT
build_classifier "338f-806r" "ACTCCTACGGGAGGCAGCAG" "GGACTACHVGGGTWTCTAAT"

# V1-V2 (27F - 338R)
# 常用序列: 27F: AGAGTTTGATCCTGGCTCAG / 338R: TGCTGCCTCCCGTAGGAGT
build_classifier "27f-338r" "AGAGTTTGATCCTGGCTCAG" "TGCTGCCTCCCGTAGGAGT"

# V4 (515F - 805R)
# 常用序列: 515F: GTGCCAGCMGCCGCGGTAA / 805R: GGACTACHVGGGTWTCTAAT
build_classifier "515f-805r" "GTGCCAGCMGCCGCGGTAA" "GGACTACHVGGGTWTCTAAT"

echo "========================================================"
echo "所有资源准备完成！分类器列表："
ls -lh *classifier.qza
echo "========================================================"