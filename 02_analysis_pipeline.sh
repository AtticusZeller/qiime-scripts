#!/bin/bash

# =======================================================
# 脚本名称: 02_analysis_pipeline.sh
# 功能: 16S 分析全流程 (自适应 单端/双端 & 多种引物区域)
# 用法: bash 02_analysis_pipeline.sh <path/to/id_download.list>
# =======================================================

set -e # 遇到错误停止

# --- 1. 参数接收与基础配置 ---
INPUT_LIST="$1"

if [ -z "$INPUT_LIST" ]; then
    echo "[ERROR] 请指定 input list 文件!"
    exit 1
fi

# 获取文件名和 ID
FILENAME=$(basename "$INPUT_LIST")
if [[ "$FILENAME" == *"_download.list" ]]; then
    PROJECT_ID="${FILENAME%_download.list}"
else
    PROJECT_ID="${FILENAME%.*}"
fi

# 工作目录
WORK_DIR="$PWD/${PROJECT_ID}_output"
mkdir -p "$WORK_DIR"
DB_BASE_DIR="$PWD/database" # 请确保指向正确的数据库目录

# 线程数
THREADS=8 
START_STEP=1 
FORCE_RUN=false 

# --- 2. 自动检测 测序类型 (Single/Dual) 和 引物区域 ---

# 2.1 检测单双端 (根据路径中是否包含 single 或 dual，如果没有则默认 dual)
if [[ "$INPUT_LIST" == *"single"* ]]; then
    DATA_LAYOUT="Single"
    echo "[CONFIG] 检测到测序类型: 单端 (Single-End)"
else
    DATA_LAYOUT="Paired"
    echo "[CONFIG] 检测到测序类型: 双端 (Paired-End)"
fi

# 2.2 检测引物区域并定义引物序列 (用于 Cutadapt)
# 注意：这里的引物序列要和 01 脚本里的一致
if [[ "$INPUT_LIST" == *"v1-v2"* ]]; then
    PRIMER_TYPE="27f-338r"
    PRIMER_F="AGAGTTTGATCCTGGCTCAG"
    PRIMER_R="TGCTGCCTCCCGTAGGAGT"
elif [[ "$INPUT_LIST" == *"v3-v4"* ]]; then
    PRIMER_TYPE="338f-806r"
    PRIMER_F="ACTCCTACGGGAGGCAGCAG"
    PRIMER_R="GGACTACHVGGGTWTCTAAT"
elif [[ "$INPUT_LIST" == *"v4"* ]]; then
    PRIMER_TYPE="515f-805r"
    PRIMER_F="GTGCCAGCMGCCGCGGTAA"
    PRIMER_R="GGACTACHVGGGTWTCTAAT"
else
    PRIMER_TYPE="338f-806r"
    PRIMER_F="ACTCCTACGGGAGGCAGCAG"
    PRIMER_R="GGACTACHVGGGTWTCTAAT"
    echo "[WARNING] 无法识别区域，默认使用 V3-V4 ($PRIMER_TYPE)"
fi

CLASSIFIER_PATH="$DB_BASE_DIR/silva-138.1-ssu-nr99-${PRIMER_TYPE}-classifier.qza"
echo "[CONFIG] 引物区域: $PRIMER_TYPE"
echo "[CONFIG] 分类器路径: $CLASSIFIER_PATH"

# --- 环境激活 ---
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.5 || { echo "环境激活失败"; exit 1; }
export PATH=~/.aspera/connect/bin:$PATH

should_run() {
    local step_num=$1
    local output_file=$2
    local step_desc=$3
    if [ "$step_num" -lt "$START_STEP" ]; then return 1; fi
    if [ -f "$output_file" ] && [ "$FORCE_RUN" = "false" ]; then
        echo "[SKIP] [$PROJECT_ID] 步骤 $step_num: $step_desc (已完成)"
        return 1
    fi
    echo "[RUN] >>> [$PROJECT_ID] 步骤 $step_num: $step_desc <<<"
    return 0
}

# ================= 流程开始 =================

# --- 步骤 1: 数据下载 (通用) ---
STEP_1_MARKER="$WORK_DIR/fastq_download_done.flag"
if should_run 1 "$STEP_1_MARKER" "数据下载"; then
    mkdir -p "$WORK_DIR/fastq"
    ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 300M -T -P33001 -k1 \
         --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp \
         --file-list "$INPUT_LIST" "$WORK_DIR/fastq"
    touch "$STEP_1_MARKER"
fi

# --- 步骤 2: 导入 QIIME2 (区分单双端) ---
DEMUX_QZA="$WORK_DIR/demux.qza"
if should_run 2 "$DEMUX_QZA" "导入数据 ($DATA_LAYOUT)"; then
    cd "$WORK_DIR"
    rm -rf new_fq && mkdir -p new_fq
    
    if [ "$DATA_LAYOUT" == "Paired" ]; then
        # === 双端导入逻辑 ===
        # 1. 整理文件名 (保留 _1/_2)
        for f in fastq/*_1.fastq.gz; do
            [ -e "$f" ] || continue
            id=$(basename "$f" | sed 's/_1\.fastq\.gz//')
            cp "$f" "new_fq/${id}_R1.fastq.gz"
        done
        for f in fastq/*_2.fastq.gz; do
            [ -e "$f" ] || continue
            id=$(basename "$f" | sed 's/_2\.fastq\.gz//')
            cp "$f" "new_fq/${id}_R2.fastq.gz"
        done
        
        # 2. 生成双端 Manifest
        ls new_fq/*_R1.fastq.gz | sed 's/_R1\.fastq\.gz//' | xargs -n 1 basename | \
        awk -v pwd="$PWD/new_fq" 'BEGIN{OFS="\t";print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"}
        { print $1, pwd"/"$1"_R1.fastq.gz", pwd"/"$1"_R2.fastq.gz" }' > manifest
        
        # 3. 导入
        IMPORT_TYPE="SampleData[PairedEndSequencesWithQuality]"
        IMPORT_FMT="PairedEndFastqManifestPhred33V2"

    else
        # === 单端导入逻辑 ===
        # 1. 整理文件名 (所有 .fastq.gz 都视为单端文件)
        # 注意: 这里假设下载的文件名直接可以作为 ID，或者通过去掉后缀作为 ID
        # 如果单端文件也带 _1 (例如 fastq-dump 出来的)，下面的逻辑也能处理
        for f in fastq/*.fastq.gz; do
            [ -e "$f" ] || continue
            # 尝试去掉 _1.fastq.gz (如果有) 或者只去掉 .fastq.gz
            filename=$(basename "$f")
            if [[ "$filename" == *"_1.fastq.gz" ]]; then
                id="${filename%_1.fastq.gz}"
            else
                id="${filename%.fastq.gz}"
            fi
            cp "$f" "new_fq/${id}.fastq.gz"
        done
        
        # 2. 生成单端 Manifest
        ls new_fq/*.fastq.gz | sed 's/\.fastq\.gz//' | xargs -n 1 basename | \
        awk -v pwd="$PWD/new_fq" 'BEGIN{OFS="\t";print "sample-id\tabsolute-filepath"}
        { print $1, pwd"/"$1".fastq.gz" }' > manifest
        
        # 3. 导入
        IMPORT_TYPE="SampleData[SequencesWithQuality]"
        IMPORT_FMT="SingleEndFastqManifestPhred33V2"
    fi

    # 检查 new_fq 是否有文件
    if [ -z "$(ls -A new_fq)" ]; then echo "[ERROR] new_fq 为空，文件整理失败"; exit 1; fi

    echo "[LOG] Importing as $IMPORT_TYPE..."
    qiime tools import \
        --type "$IMPORT_TYPE" \
        --input-path manifest \
        --output-path demux.qza \
        --input-format "$IMPORT_FMT"

    qiime demux summarize --i-data demux.qza --o-visualization demux-summary.qzv
    cd ..
fi

# --- 步骤 3: Cutadapt 去引物 (区分单双端) ---
TRIMMED_QZA="$WORK_DIR/primer-trimmed-demux.qza"
if should_run 3 "$TRIMMED_QZA" "Cutadapt 去引物"; then
    if [ "$DATA_LAYOUT" == "Paired" ]; then
        # 双端切引物
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences "$DEMUX_QZA" \
            --p-cores $THREADS \
            --p-no-indels \
            --p-front-f "$PRIMER_F" \
            --p-front-r "$PRIMER_R" \
            --o-trimmed-sequences "$TRIMMED_QZA"
    else
        # 单端切引物 (只切 5' 端，即 Front)
        qiime cutadapt trim-single \
            --i-demultiplexed-sequences "$DEMUX_QZA" \
            --p-cores $THREADS \
            --p-no-indels \
            --p-front "$PRIMER_F" \
            --o-trimmed-sequences "$TRIMMED_QZA"
    fi
    # 可视化
    qiime demux summarize --i-data "$TRIMMED_QZA" --o-visualization "$WORK_DIR/primer-trimmed-demux.qzv"
fi

# --- 步骤 4: 双端合并 (单端跳过) ---
# 定义 PRE_FILTER_QZA 变量，指向进入质控步骤的文件
if [ "$DATA_LAYOUT" == "Paired" ]; then
    JOINED_QZA="$WORK_DIR/demux-joined.qza"
    if should_run 4 "$JOINED_QZA" "双端合并"; then
        qiime vsearch merge-pairs \
            --i-demultiplexed-seqs "$TRIMMED_QZA" \
            --o-unmerged-sequences "$WORK_DIR/unmerged_demux-joined.qza" \
            --p-threads $THREADS \
            --o-merged-sequences "$JOINED_QZA"
        qiime demux summarize --i-data "$JOINED_QZA" --o-visualization "$WORK_DIR/demux-joined-summary.qzv"
    fi
    PRE_FILTER_QZA="$JOINED_QZA"
else
    # 单端模式：直接跳过合并，Step 3 的产物直接进入 Step 5
    echo "[INFO] 单端数据，跳过双端合并步骤。"
    PRE_FILTER_QZA="$TRIMMED_QZA"
fi

# --- 步骤 5: 质控过滤 (通用) ---
FILTERED_QZA="$WORK_DIR/demux-filtered.qza"
if should_run 5 "$FILTERED_QZA" "质控过滤"; then
    # 注意：这里的输入是 PRE_FILTER_QZA (可能是 joined 也可能是 trimmed-single)
    qiime quality-filter q-score \
        --p-min-quality 20 \
        --i-demux "$PRE_FILTER_QZA" \
        --o-filtered-sequences "$FILTERED_QZA" \
        --o-filter-stats "$WORK_DIR/filter-stats.qza"
fi

# --- 步骤 6: Deblur 去噪 (自动计算长度) ---
TABLE="$WORK_DIR/feature-table.qza"
REP_SEQS="$WORK_DIR/repset-seqs.qza"

if should_run 6 "$TABLE" "Deblur 去噪"; then
    echo "[AUTO] 正在计算 Deblur trim-length..."
    TEMP_EXPORT_DIR="$WORK_DIR/temp_len_check"
    rm -rf "$TEMP_EXPORT_DIR"
    qiime tools export --input-path "$FILTERED_QZA" --output-path "$TEMP_EXPORT_DIR"
    
    # 查找 fastq 文件
    SAMPLE_FQ=$(find "$TEMP_EXPORT_DIR" -name "*.fastq.gz" | head -n 1)
    
    # 计算第 10% 分位数的长度
    TRIM_LEN=$(zcat "$SAMPLE_FQ" | head -n 16000 | awk 'NR%4==2 {print length($0)}' | sort -n | awk 'BEGIN{c=0} {a[c]=$1; c++} END{print a[int(c*0.1)]}')
    
    # 安全检查
    if [ -z "$TRIM_LEN" ] || [ "$TRIM_LEN" -lt 100 ]; then
        echo "[WARNING] 计算长度异常 ($TRIM_LEN)，重置为 200"
        TRIM_LEN=200
    fi
    echo "[AUTO] 截断长度: ${TRIM_LEN} bp"
    
    qiime deblur denoise-16S \
        --i-demultiplexed-seqs "$FILTERED_QZA" \
        --p-trim-length $TRIM_LEN \
        --p-sample-stats \
        --p-jobs-to-start $THREADS \
        --p-min-reads 1 \
        --o-representative-sequences "$REP_SEQS" \
        --o-table "$TABLE" \
        --o-stats "$WORK_DIR/deblur-stats.qza"
        
    rm -rf "$TEMP_EXPORT_DIR"
fi

# --- 步骤 7, 8, 9 (通用) ---
# 后续分类、折叠、导出逻辑对于单双端完全一致

TAXONOMY="$WORK_DIR/taxonomy.qza"
if should_run 7 "$TAXONOMY" "物种注释"; then
    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER_PATH" \
        --i-reads "$REP_SEQS" \
        --p-n-jobs $THREADS \
        --o-classification "$TAXONOMY"
fi

GENUS_TABLE="$WORK_DIR/table-L6-genus.qza"
if should_run 8 "$GENUS_TABLE" "Collapse Genus"; then
    # # === 新增：自动检查分类结果是否正常 ===
    # echo "[CHECK] 正在检查分类结果..."
    # qiime metadata tabulate \
    #     --m-input-file "$TAXONOMY" \
    #     --o-visualization "$WORK_DIR/taxonomy_check.qzv"
    
    # # 导出 taxonomy 文本查看前几行
    # qiime tools export --input-path "$TAXONOMY" --output-path "$WORK_DIR/temp_tax_check"
    
    # echo "--- 分类结果预览 (前 5 行) ---"
    # head -n 5 "$WORK_DIR/temp_tax_check/taxonomy.tsv"
    # echo "------------------------------"
    
    # # 检查是否包含 Unassigned
    # UNASSIGNED_COUNT=$(grep -c "Unassigned" "$WORK_DIR/temp_tax_check/taxonomy.tsv" || true)
    # TOTAL_COUNT=$(wc -l < "$WORK_DIR/temp_tax_check/taxonomy.tsv")
    
    # if [ "$UNASSIGNED_COUNT" -gt $((TOTAL_COUNT / 2)) ]; then
    #     echo "[WARNING] 警告：超过 50% 的序列未被分类 (Unassigned)。"
    #     echo "可能原因：1. 引物区域不匹配 (已更换全长分类器解决) 2. 序列方向反了 3. 数据非 16S"
    # fi
    qiime taxa collapse \
        --i-table "$TABLE" \
        --i-taxonomy "$TAXONOMY" \
        --p-level 6 \
        --o-collapsed-table "$GENUS_TABLE"
fi

FINAL_TSV="$WORK_DIR/${PROJECT_ID}_genus-table.tsv"
if should_run 9 "$FINAL_TSV" "导出 TSV"; then
    qiime tools export --input-path "$GENUS_TABLE" --output-path "$WORK_DIR/exported-genus-table"
    biom convert -i "$WORK_DIR/exported-genus-table/feature-table.biom" -o "$FINAL_TSV" --to-tsv
    echo "[SUCCESS] [$PROJECT_ID] 处理完成！TSV: $FINAL_TSV"
fi