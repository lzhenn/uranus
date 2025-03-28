import os
import time
from uranus.lib import utils
from datetime import datetime, timedelta
# Print the contents of the precfg file
def print_cfg(precfg):
    for section in precfg.sections():
        print(f"[{section}]")
        for option in precfg.options(section):
            print(f"{option} = {precfg.get(section, option)}")

# Function to read the last line of the file
def read_last_line(file_path):
    with open(file_path, 'rb') as f:
        f.seek(-2, os.SEEK_END)  # Move to the end of the file
        while f.read(1) != b'\n':  # Go backwards until we find a newline
            f.seek(-2, os.SEEK_CUR)
        return f.readline().decode()

def monitor_file(rsl_dir, trigger_time):
    
    # Loop until the trigger time is reached 
    rsl_file=os.path.join(rsl_dir,'rsl.out.0000')
    last_line = None
    unchanged_count = 0

    while True:
        if not(os.path.exists(rsl_file)):
            utils.write_log("rsl.out.0000 does not exist. waiting for 30 seconds...")
            time.sleep(30)
            continue
        current_line = read_last_line(rsl_file)
        
        
        curr_time = datetime.now()  # Example trigger time, set as needed
        utils.write_log(f"Current Line: {current_line}")
        if current_line != last_line:
            last_line = current_line
            unchanged_count = 0  # Reset count if line changes
        else:
            unchanged_count += 1

        # Check for unchanged line for 5 minutes
        if unchanged_count >= 30:  # 30 checks of 10 seconds each = 5 minutes
            utils.write_log("The last line has not changed for 5 minutes. WRF.exe may abort.")
            return -1
        
        # Check for the specific timestamp
        if "Timing for main: time" in current_line:
            timestamp_str = current_line.split("time ")[1].strip()
            timestamp_str = timestamp_str.split(" ")[0].strip()
            timestamp = datetime.strptime(timestamp_str, '%Y-%m-%d_%H:%M:%S')
            utils.write_log(f">>>wrf.exe<<< Running @ {timestamp}")
            # Compare with trigger time
            if abs((timestamp - trigger_time).total_seconds()) <= 3000:  # 50 minutes
                return 0 

        time.sleep(5)  # Wait for 10 seconds before the next check

