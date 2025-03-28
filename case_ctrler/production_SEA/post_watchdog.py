import os
import time
import re
from multiprocessing import Pool
from pathlib import Path
import subprocess
import logging
from datetime import datetime

# 配置参数
MONITOR_DIR = "/home/lzhenn/hqnfs/WRF-4.1.5_org_P2/run"  # 监控目录
CACHE_DIR = "/home/lzhenn/hqnfs/temp/cache"          # 缓存目录
NTASKS = 4                                           # 最大并行进程数
TIMEOUT = 7200                                       # 1小时超时（秒）
WRFOUT_PATTERN = r"wrfout_d\d{2}_\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}"  # wrfout文件格式
#WRFRST_PATTERN = r"wrfrst_d\d{2}_\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}"  # wrfrst文件格式
WRFRST_PATTERN = r"wrfrst_d\d{2}_\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}_\d{4}"  # wrfrst文件格式

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 确保缓存目录存在
Path(CACHE_DIR).mkdir(parents=True, exist_ok=True)

def parse_timestamp(filename, pattern):
    """从文件名中提取时间戳并转换为datetime对象"""
    match = re.search(r"(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})", filename)
    if match:
        timestamp_str = match.group(1)  # 直接捕获 YYYY-MM-DD_HH:MM:SS
        return datetime.strptime(timestamp_str, "%Y-%m-%d_%H:%M:%S")
    logger.warning(f"Could not parse timestamp from filename: {filename}")
    return None
def convert_nc3_to_nc4(file_path):
    """将nc3文件转换为nc4并处理"""
    try:
        file_name = os.path.basename(file_path)
        output_path = os.path.join(CACHE_DIR, f"{file_name}")
        
        # 使用ncks命令转换格式
        cmd = ["ncks", "-4", "-L", "1", file_path, output_path]
        subprocess.run(cmd, check=True)
        
        # 验证转换成功后删除原文件
        if os.path.exists(output_path):
            os.remove(file_path)
            logger.info(f"Converted and removed: {file_path}")
        else:
            logger.error(f"Conversion failed for: {file_path}")
    except Exception as e:
        logger.error(f"Error converting {file_path}: {str(e)}")

def get_wrfout_files(directory):
    """获取目录下符合条件的wrfout文件列表，并按时间戳排序"""
    files = [f for f in os.listdir(directory) 
            if re.match(WRFOUT_PATTERN, f) and os.path.isfile(os.path.join(directory, f))]
    
    # 按时间戳排序
    files_with_paths = [os.path.join(directory, f) for f in files]
    files_with_paths.sort(key=lambda x: parse_timestamp(os.path.basename(x), WRFOUT_PATTERN))
    return files_with_paths

def get_wrfrst_files(directory):
    """获取目录下符合条件的wrfrst文件列表，并按时间戳排序"""
    files = [f for f in os.listdir(directory) 
            if re.match(WRFRST_PATTERN, f) and os.path.isfile(os.path.join(directory, f))]
    
    # 按时间戳排序
    files_with_paths = [os.path.join(directory, f) for f in files]
    files_with_paths.sort(key=lambda x: parse_timestamp(os.path.basename(x), WRFRST_PATTERN))
    return files_with_paths

def process_files(file_list):
    """使用多进程处理wrfout文件转换"""
    with Pool(processes=min(NTASKS, len(file_list))) as pool:
        pool.map(convert_nc3_to_nc4, file_list)

def process_wrfrst_files(wrfrst_files):
    """处理wrfrst文件：按时间戳创建目录，移动文件并压缩"""
    if not wrfrst_files:
        return
    
    # 按时间戳分组
    time_groups = {}
    for file_path in wrfrst_files:
        timestamp = parse_timestamp(os.path.basename(file_path), WRFRST_PATTERN)
        if timestamp:
            dir_name = f"wrfrst_{timestamp.strftime('%Y-%m-%d_%H')}"
            if dir_name not in time_groups:
                time_groups[dir_name] = []
            time_groups[dir_name].append(file_path)

    # 处理每个时间组
    for dir_name, files in time_groups.items():
        target_dir = os.path.join(MONITOR_DIR, dir_name)
        Path(target_dir).mkdir(exist_ok=True)
        
        # 移动文件到对应目录
        for file_path in files:
            file_name = os.path.basename(file_path)
            os.rename(file_path, os.path.join(target_dir, file_name))
            logger.info(f"Moved {file_name} to {target_dir}")
        
        # 压缩目录
        tar_name = f"{dir_name}.tar.gz"
        tar_path = os.path.join(CACHE_DIR, tar_name)
        cmd = ["tar", "-czvf", tar_path, "-C", MONITOR_DIR, dir_name]
        subprocess.run(cmd, check=True)
        logger.info(f"Compressed {dir_name} to {tar_path}")
        
        # 删除原目录
        for file in os.listdir(target_dir):
            os.remove(os.path.join(target_dir, file))
        os.rmdir(target_dir)
        logger.info(f"Removed temporary directory: {target_dir}")

def main():
    logger.info("Starting directory monitoring...")
    last_activity_time = time.time()
    
    while True:
        current_time = time.time()
        dt=current_time - last_activity_time        
        # 获取当前wrfout和wrfrst文件列表（已排序）
        wrfout_files = get_wrfout_files(MONITOR_DIR)
        wrfrst_files = get_wrfrst_files(MONITOR_DIR)
        
        # 处理wrfrst文件
        if wrfrst_files:
            logger.info(f"Found {len(wrfrst_files)} wrfrst files, wait for 600s before processing...")
            time.sleep(600)
            process_wrfrst_files(wrfrst_files)
            last_activity_time = time.time()
        
        # 处理wrfout文件
        if wrfout_files:
            if len(wrfout_files) > NTASKS:
                logger.info(f"Found {len(wrfout_files)} wrfout files (> {NTASKS}), starting processing")
                process_files(wrfout_files)
                last_activity_time = time.time()
            else:
                logger.info(f"Found {len(wrfout_files)} wrfout files (<= {NTASKS}), {dt:.0f}/{TIMEOUT} seconds elapsed...")
        # 检查是否超过1小时没有新文件
        if dt > TIMEOUT:
            logger.info("No new files for 1 hour, stopping monitor")
            break
        logger.info(f"No files to process, {dt:.0f}/{TIMEOUT} seconds elapsed...")
    
        # 每30秒检查一次
        time.sleep(30)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Monitoring stopped by user")
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
