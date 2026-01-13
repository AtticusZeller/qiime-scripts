import os
import glob
import subprocess
import concurrent.futures
import time
from pathlib import Path

# ================= 配置区域 =================
# 存放 _download.list 文件的目录 (支持 . 为当前目录)
INPUT_DIR = "./inputs"

# Shell 脚本路径
PIPELINE_SCRIPT = "./02_analysis_pipeline.sh"

# 最大并发任务数
# 警告: 每个任务内部还会调用多线程 (qiime2 threads)，请谨慎设置
# 推荐公式: 服务器总核心数 / 单个任务 internal threads (02脚本里的THREADS)
MAX_WORKERS = 8
# ===========================================

def run_task(list_file):
    """
    调用 Shell 脚本处理单个 List 文件
    """
    file_path = Path(list_file)
    file_name = file_path.name
    
    print(f"[开始] 调度任务: {file_name}")
    
    start_time = time.time()
    
    try:
        # 调用 bash 脚本
        # check=True 会在脚本返回非0状态码时抛出异常
        # capture_output=False 让脚本的输出直接打印到屏幕，方便观察（或者重定向到日志文件）
        result = subprocess.run(
            ["bash", PIPELINE_SCRIPT, str(file_path)],
            check=True,
            text=True
        )
        elapsed = time.time() - start_time
        return f"[完成] {file_name} (耗时: {elapsed:.2f}s)"
        
    except subprocess.CalledProcessError as e:
        return f"[失败] {file_name} (退出代码: {e.returncode})"
    except Exception as e:
        return f"[错误] {file_name}: {str(e)}"

def main():
    # 1. 检查输入目录
    if not os.path.exists(INPUT_DIR):
        print(f"错误: 输入目录 '{INPUT_DIR}' 不存在，请创建并放入 .list 文件")
        return

    # 2. 搜索所有 _download.list 文件
    # 假设你的文件模式是 *_download.list
    search_pattern = os.path.join(INPUT_DIR, "*_download.list")
    list_files = glob.glob(search_pattern)
    
    if not list_files:
        print(f"警告: 在 '{INPUT_DIR}' 下未找到任何匹配 '{search_pattern}' 的文件")
        return

    print(f"找到 {len(list_files)} 个任务，准备开始处理...")
    print(f"并发数: {MAX_WORKERS}")
    print("-" * 50)

    # 3. 多进程并发执行
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # 提交所有任务
        future_to_file = {executor.submit(run_task, f): f for f in list_files}
        
        # 处理结果 (as_completed 会在任务完成时立即返回)
        for future in concurrent.futures.as_completed(future_to_file):
            original_file = future_to_file[future]
            try:
                status_msg = future.result()
                print(status_msg)
            except Exception as exc:
                print(f"[异常] 文件 {original_file} 执行时发生未知错误: {exc}")

    print("-" * 50)
    print("所有任务处理完毕。")

if __name__ == "__main__":
    main()