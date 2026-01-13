#!/bin/bash

# =======================================================
# 脚本名称: 02_analysis_pipeline.sh
# 功能: 16S 分析全流程 (单样本处理版)
# 用法: bash 02_analysis_pipeline.sh <path/to/id_download.list>
# =======================================================

set -e # 遇到错误停止

# --- 1. 参数接收与 ID 解析 ---
INPUT_LIST="$1"

if [ -z "$INPUT_LIST" ]; then
    echo "[ERROR] 请指定 input list 文件!"
    echo "用法: bash $0 <list_file>"
    exit 1
fi

if [ ! -f "$INPUT_LIST" ]; then
    echo "[ERROR] 文件不存在: $INPUT_LIST"
    exit 1
fi

# 获取文件名 (不带路径)
FILENAME=$(basename "$INPUT_LIST")
# 提取 ID (假设文件名格式为 ID_download.list，如果不是，则直接用文件名做ID)
if [[ "$FILENAME" == *"_download.list" ]]; then
    PROJECT_ID="${FILENAME%_download.list}"
else
    PROJECT_ID="${FILENAME%.*}"
fi

echo "[INIT] 检测到 Project ID: $PROJECT_ID"
echo "[INIT] 输入文件: $INPUT_LIST"

# --- 2. 配置参数 ---

# 所有的输出都会在这个文件夹下
WORK_DIR="$PWD/${PROJECT_ID}_output"
mkdir -p "$WORK_DIR"

# 数据库路径 (请确保这里指向你 01 脚本生成的真实路径)
# 建议使用绝对路径，防止目录切换导致找不到
DB_BASE_DIR="$PWD/database" 
CLASSIFIER_PATH="$DB_BASE_DIR/silva-138.1-ssu-nr99-338f-806r-classifier.qza"

# 单个任务内部使用的线程数
# 注意：如果你由 03 脚本并发调用，这里不要设置太大，否则服务器会卡死
# 比如：服务器共64核，你打算同时跑4个任务，那这里设为 16
THREADS=8

# --- 控制流程 ---
START_STEP=5
FORCE_RUN=true

# --- 环境激活 ---
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.5 || { echo "环境激活失败"; exit 1; }
export PATH=~/.aspera/connect/bin:$PATH

# 辅助函数
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

# --- 步骤 1: 数据下载 ---
STEP_1_MARKER="$WORK_DIR/fastq_download_done.flag"
if should_run 1 "$STEP_1_MARKER" "数据下载"; then
    mkdir -p "$WORK_DIR/fastq"
    # 注意：这里使用传入的 INPUT_LIST
    ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 300M -T -P33001 -k1 \
         --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp \
         --file-list "$INPUT_LIST" "$WORK_DIR/fastq"
    touch "$STEP_1_MARKER"
fi

# --- 步骤 2: 导入 QIIME2 ---
DEMUX_QZA="$WORK_DIR/demux.qza"
if should_run 2 "$DEMUX_QZA" "导入数据"; then
    cd "$WORK_DIR"
    
    # 2.1 清理并准备目录
    if [ -d "new_fq" ]; then rm -rf "new_fq"; fi
    mkdir -p new_fq
    
    echo "[LOG] 正在整理 Fastq 文件..."
    
    # 2.2 更加稳健的文件重命名与移动逻辑
    # 仅处理以 _1.fastq.gz 或 _2.fastq.gz 结尾的文件，自动忽略不合规文件
    
    # 处理 R1 (Forward)
    for f in fastq/*_1.fastq.gz; do
        [ -e "$f" ] || continue # 如果没有匹配文件则跳过
        filename=$(basename "$f")
        # 提取 ID: 去掉后缀 _1.fastq.gz
        id="${filename%_1.fastq.gz}" 
        cp "$f" "new_fq/${id}_R1.fastq.gz"
    done

    # 处理 R2 (Reverse)
    for f in fastq/*_2.fastq.gz; do
        [ -e "$f" ] || continue
        filename=$(basename "$f")
        # 提取 ID: 去掉后缀 _2.fastq.gz
        id="${filename%_2.fastq.gz}"
        cp "$f" "new_fq/${id}_R2.fastq.gz"
    done

    # 检查 new_fq 是否有文件
    count=$(ls new_fq/*.fastq.gz 2>/dev/null | wc -l)
    if [ "$count" -eq 0 ]; then
        echo "[ERROR] new_fq 目录为空！未找到 _1.fastq.gz / _2.fastq.gz 结尾的文件。"
        echo "[HINT] 请检查 fastq 文件夹下的文件名格式。"
        exit 1
    fi

    # 2.3 生成 Manifest (基于整理后的 new_fq 目录)
    echo "[LOG] 生成 Manifest 文件..."
    
    # 逻辑：
    # 1. 列出 new_fq 下所有文件
    # 2. 移除 _R1.fastq.gz 或 _R2.fastq.gz 得到纯 SampleID
    # 3. 去重
    # 4. 使用 awk 拼接绝对路径
    ls new_fq/*.fastq.gz | sed -E 's/_R[12]\.fastq\.gz//g' | xargs -n 1 basename | uniq | \
    awk -v pwd="$PWD/new_fq" 'BEGIN{OFS="\t";print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"}
    {
        # 构造预期的 R1 和 R2 路径
        r1 = pwd "/" $1 "_R1.fastq.gz"
        r2 = pwd "/" $1 "_R2.fastq.gz"
        
        # 这一步是用来生成 manifest 内容的，QIIME2 导入时会检查文件是否存在
        print $1, r1, r2
    }' > manifest

    echo "[LOG] 正在导入 QIIME 2..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path manifest \
        --output-path demux.qza \
        --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize --i-data demux.qza --o-visualization demux-summary.qzv
    
    cd .. # 回到上级目录
fi
# --- 步骤 3: Cutadapt ---
TRIMMED_QZA="$WORK_DIR/primer-trimmed-demux.qza"
if should_run 3 "$TRIMMED_QZA" "Cutadapt 去引物"; then
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "$DEMUX_QZA" \
        --p-cores $THREADS \
        --p-no-indels \
        --p-front-f ACTCCTACGGGAGGCAGCAG \
        --p-front-r GGACTACHVGGGTWTCTAAT \
        --o-trimmed-sequences "$TRIMMED_QZA" \
        --verbose > "$WORK_DIR/cutadapt_log.txt" 2>&1
fi

# --- 步骤 4: 双端合并 ---
JOINED_QZA="$WORK_DIR/demux-joined.qza"
if should_run 4 "$JOINED_QZA" "双端合并"; then
    qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "$TRIMMED_QZA" \
        --o-unmerged-sequences "$WORK_DIR/unmerged_demux-joined.qza" \
        --p-threads $THREADS \
        --o-merged-sequences "$JOINED_QZA"
    
    qiime demux summarize --i-data "$JOINED_QZA" --o-visualization "$WORK_DIR/demux-joined-summary.qzv"
fi

# --- 步骤 5: 质控过滤 ---
FILTERED_QZA="$WORK_DIR/demux-joined-filtered.qza"
if should_run 5 "$FILTERED_QZA" "质控过滤"; then
    qiime quality-filter q-score \
        --p-min-quality 20 \
        --i-demux "$JOINED_QZA" \
        --o-filtered-sequences "$FILTERED_QZA" \
        --o-filter-stats "$WORK_DIR/demux-joined-filter-stats.qza"
fi

# --- 步骤 6: Deblur 去噪 (自动计算 trim-length) ---
TABLE="$WORK_DIR/feature-table.qza"
REP_SEQS="$WORK_DIR/repset-seqs.qza"

if should_run 6 "$TABLE" "Deblur 去噪 (自动计算长度)"; then
    echo "[AUTO] 正在根据数据分布自动计算 Deblur trim-length..."
    
    # 1. 临时导出过滤后的数据以检测长度
    TEMP_EXPORT_DIR="$WORK_DIR/temp_len_check"
    rm -rf "$TEMP_EXPORT_DIR" # 清理旧的
    qiime tools export \
        --input-path "$FILTERED_QZA" \
        --output-path "$TEMP_EXPORT_DIR"
    
    # 2. 找到第一个样本的 fastq 文件
    SAMPLE_FQ=$(find "$TEMP_EXPORT_DIR" -name "*.fastq.gz" | head -n 1)
    
    if [ -z "$SAMPLE_FQ" ]; then
        echo "[ERROR] 无法找到导出的 Fastq 文件，请检查上一步质控结果是否为空！"
        exit 1
    fi
    
    # 3. 计算 trim-length
    # 逻辑：取前 4000 条序列 -> 计算长度 -> 排序 -> 取第 10% 分位数的长度
    # 这意味着 90% 的序列长度都大于等于这个值，这样可以保留绝大多数数据
    TRIM_LEN=$(zcat "$SAMPLE_FQ" | head -n 16000 | awk 'NR%4==2 {print length($0)}' | sort -n | awk 'BEGIN{c=0} {a[c]=$1; c++} END{print a[int(c*0.1)]}')
    
    # 防止计算出异常值 (例如太短)
    if [ "$TRIM_LEN" -lt 100 ]; then
        echo "[WARNING] 计算出的截断长度过短 ($TRIM_LEN)，重置为 250 (V4保守值) 或请检查数据。"
        TRIM_LEN=250
    fi
    
    echo "[AUTO] 自动确定的截断长度为: ${TRIM_LEN} bp (基于第 10 百分位)"

    # 4. 运行 Deblur
    qiime deblur denoise-16S \
        --i-demultiplexed-seqs "$FILTERED_QZA" \
        --p-trim-length $TRIM_LEN \
        --p-sample-stats \
        --p-jobs-to-start $THREADS \
        --p-min-reads 1 \
        --o-representative-sequences "$REP_SEQS" \
        --o-table "$TABLE" \
        --o-stats "$WORK_DIR/deblur-stats.qza"
    
    # 清理临时文件
    rm -rf "$TEMP_EXPORT_DIR"
    
    # 可视化
    qiime feature-table summarize --i-table "$TABLE" --o-visualization "$WORK_DIR/feature-table.qzv"
    qiime deblur visualize-stats --i-deblur-stats "$WORK_DIR/deblur-stats.qza" --o-visualization "$WORK_DIR/deblur-stats.qzv"
fi

# --- 步骤 7: 物种注释 ---
TAXONOMY="$WORK_DIR/taxonomy.qza"
if should_run 7 "$TAXONOMY" "物种注释"; then
    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER_PATH" \
        --i-reads "$REP_SEQS" \
        --p-n-jobs $THREADS \
        --o-classification "$TAXONOMY"
fi

# --- 步骤 8: Collapse Genus ---
GENUS_TABLE="$WORK_DIR/table-L6-genus.qza"
if should_run 8 "$GENUS_TABLE" "Collapse Genus"; then
    # === 新增：自动检查分类结果是否正常 ===
    echo "[CHECK] 正在检查分类结果..."
    qiime metadata tabulate \
        --m-input-file "$TAXONOMY" \
        --o-visualization "$WORK_DIR/taxonomy_check.qzv"
    
    # 导出 taxonomy 文本查看前几行
    qiime tools export --input-path "$TAXONOMY" --output-path "$WORK_DIR/temp_tax_check"
    
    echo "--- 分类结果预览 (前 5 行) ---"
    head -n 5 "$WORK_DIR/temp_tax_check/taxonomy.tsv"
    echo "------------------------------"
    
    # 检查是否包含 Unassigned
    UNASSIGNED_COUNT=$(grep -c "Unassigned" "$WORK_DIR/temp_tax_check/taxonomy.tsv" || true)
    TOTAL_COUNT=$(wc -l < "$WORK_DIR/temp_tax_check/taxonomy.tsv")
    
    if [ "$UNASSIGNED_COUNT" -gt $((TOTAL_COUNT / 2)) ]; then
        echo "[WARNING] 警告：超过 50% 的序列未被分类 (Unassigned)。"
        echo "可能原因：1. 引物区域不匹配 (已更换全长分类器解决) 2. 序列方向反了 3. 数据非 16S"
    fi
    qiime taxa collapse \
        --i-table "$TABLE" \
        --i-taxonomy "$TAXONOMY" \
        --p-level 6 \
        --o-collapsed-table "$GENUS_TABLE"
fi

# --- 步骤 9: 导出 TSV ---
# 这里的输出文件名包含了 ID
FINAL_TSV="$WORK_DIR/${PROJECT_ID}_genus-table.tsv"
if should_run 9 "$FINAL_TSV" "导出 TSV"; then
    qiime tools export \
        --input-path "$GENUS_TABLE" \
        --output-path "$WORK_DIR/exported-genus-table"
    
    biom convert \
        -i "$WORK_DIR/exported-genus-table/feature-table.biom" \
        -o "$FINAL_TSV" \
        --to-tsv
    
    echo "[SUCCESS] [$PROJECT_ID] 处理完成！结果: $FINAL_TSV"
fi