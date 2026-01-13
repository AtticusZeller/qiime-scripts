import csv
import sys
import re
import glob
from pathlib import Path

def process_tsv(input_file):
    # 1. 确定输出文件名 (正则提取 PRJ 号，或者使用文件名)
    match = re.search(r'(PRJ[A-Z0-9]+)', input_file.name)
    prefix = match.group(1) if match else input_file.stem
    output_file = f"{prefix}_download.list"

    links = []
    target_col = 'fastq_aspera' # 这里可以修改目标列

    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if target_col not in reader.fieldnames:
                print(f"[跳过] {input_file}: 找不到列 {target_col}")
                return

            for row in reader:
                if row[target_col]:
                    # 分割分号，并去除冒号前缀 (split(':')[-1])
                    raw_paths = row[target_col].split(';')
                    clean_paths = [p.split(':')[-1].strip() for p in raw_paths if p.strip()]
                    links.extend(clean_paths)

        if links:
            with open(output_file, 'w', encoding='utf-8') as out:
                out.write('\n'.join(links) + '\n')
            print(f"[成功] {input_file} -> {output_file} ({len(links)} 条链接)")
        else:
            print(f"[警告] {input_file} 内容为空")

    except Exception as e:
        print(f"[错误] 处理 {input_file} 失败: {e}")

def main():
    # 如果有命令行参数则使用参数，否则扫描当前目录 *.tsv
    files = sys.argv[1:] if len(sys.argv) > 1 else glob.glob("*.tsv")

    if not files:
        print("未找到 TSV 文件")
        return

    for f in files:
        process_tsv(Path(f))

if __name__ == "__main__":
    main()