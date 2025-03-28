import subprocess
import time
from threading import Thread

def check_process(process):
    # 等待B.py进程完成
    process.wait()
    print("runtime completed")

def main():
    # 使用subprocess启动B.py
    process = subprocess.Popen(['python', 'wrf_runtime.py'])
    
    # 创建线程来监控B.py的完成状态
    monitor_thread = Thread(target=check_process, args=(process,))
    monitor_thread.start()
    
    # A.py继续执行其他任务
    print("A.py continuing execution...")
    for i in range(5):
        print(f"A.py working: {i}")
        time.sleep(1)
    
    # 主线程可以继续其他工作
    # monitor_thread会在B.py完成后打印"B completed"

if __name__ == "__main__":
    main()