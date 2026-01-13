#!/bin/bash

# =======================================================
# 脚本名称: 04_create_visualization.sh
# 功能: 生成可交互的物种堆叠图 (Taxa Barplot)
# 用法: bash 04_create_visualization.sh <output_directory>
# 示例: bash 04_create_visualization.sh PRJNA770746_output
# =======================================================

set -e

WORK_DIR="$1"

# 1. 检查参数
if [ -z "$WORK_DIR" ]; then
    echo "错误: 请指定输出目录路径"
    echo "用法: bash $0 <output_folder>"
    exit 1
fi

# 2. 检查必要文件是否存在
# 我们需要: 
# 1. 原始特征表 (feature-table.qza) - 注意是未经 collapse 的
# 2. 物种注释 (taxonomy.qza)
# 3. 样本元数据 (manifest) - barplot 需要元数据来排列样本

TABLE_QZA="$WORK_DIR/feature-table.qza"
TAXONOMY_QZA="$WORK_DIR/taxonomy.qza"
MANIFEST="$WORK_DIR/manifest"

if [ ! -f "$TABLE_QZA" ] || [ ! -f "$TAXONOMY_QZA" ]; then
    echo "错误: 在目录 $WORK_DIR 中找不到 feature-table.qza 或 taxonomy.qza"
    echo "请确保 02 脚本已运行到步骤 7。"
    exit 1
fi

echo "[LOG] 正在生成交互式物种堆叠图..."

# # 3. 激活环境 (如果尚未激活)
# source ~/.bashrc
# eval "$(conda shell.bash hook)"
# conda activate qiime2-amplicon-2024.5

# 4. 生成可视化文件 (.qzv)
# --m-metadata-file 使用之前的 manifest 文件，因为它包含了 sample-id
qiime taxa barplot \
    --i-table "$TABLE_QZA" \
    --i-taxonomy "$TAXONOMY_QZA" \
    --m-metadata-file "$MANIFEST" \
    --o-visualization "$WORK_DIR/taxa-bar-plots.qzv"

echo "========================================================"
echo "[SUCCESS] 可视化文件生成成功!"
echo "文件路径: $WORK_DIR/taxa-bar-plots.qzv"
echo "--------------------------------------------------------"
echo "操作指南:"
echo "1. 下载该 .qzv 文件到本地电脑。"
echo "2. 打开网站: https://view.qiime2.org"
echo "3. 将文件拖入网页。"
echo "4. 在网页顶部 'Taxonomic Level' 下拉菜单选择 'Level 7 (Species)'。"
echo "5. 如果你想导出表格:"
echo "   - 点击顶部的 'CSV' 按钮下载当前层级的表格。"
echo "========================================================"