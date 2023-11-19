#/usr/bin/env python3
"""
    Commonly used utilities

    Function    
    ---------------
    throw_error(msg):
        throw error and exit
    
    write_log(msg, lvl=20):
        write logging log to log file
    
"""
# ---imports---
import logging
import pkg_resources
from . import const
import shutil
import datetime
from tempfile import mkstemp

def build_wrfcmd(mach_name, bashrc, exeroot, mpicmd, ntasks, exename):
    mach_dic=const.MACHINE_DIC[mach_name]
    nnodes=max(1,ntasks//mach_dic['corespernode'])
    if 'hqlx' in mach_name:
        return f'ssh {mach_name} "source {bashrc}; cd {exeroot};{mpicmd} -np {ntasks} ./{exename}"'
    elif mach_name == 'th2':
        return f'source {bashrc}; cd {exeroot};{mpicmd} -N {nnodes} -n {ntasks} ./{exename}'
    elif mach_name == 'pird':
        return f'source {bashrc}; cd {exeroot};{mpicmd} -N {nnodes} -n {ntasks} ./{exename}.sh'
    else:
        return f'source {bashrc}; cd {exeroot};{mpicmd} -np {ntasks} ./{exename}'
    
def parse_init_time(tstr):
    try:
        ts=datetime.datetime.strptime(
            tstr,'%Y%m%d%H')
    except:
        if tstr.startswith('T-'):
            num_days = int(tstr[2:])
            today = datetime.datetime.today()
            today = today.replace(hour=0,minute=0,second=0,microsecond=0)
            ts=today - datetime.timedelta(days=num_days)
    return ts
def parse_fmt_timepath(tgt_time, fmtpath):
    '''
    parse time string to datetime object
    '''
    seg_path=fmtpath.split('@')
    parsed_path=''
    for seg in seg_path:
        if seg.startswith('%'):
            parsed_path+=tgt_time.strftime(seg)
        else:
            parsed_path+=seg
    return parsed_path

def sed_wrf_timeline(pattern, ts, source, fmt="'%Y-%m-%d_%H:%M:%S'", dedest=None):
    ts_list=[ts.strftime(fmt)]*4
    temp_str=','.join(ts_list)
    sedline(pattern,f'{pattern} = {temp_str}',source)
def sedline(pattern, replace, source, dest=None, count=1):
    """Reads a source file and writes the destination file.

    In each line, replaces pattern with replace.

    Args:
        pattern (str): pattern to match (can be re.pattern)
        replace (str): replacement str
        source  (str): input filename
        count (int): number of occurrences to replace
        dest (str):   destination filename, if not given, source will be over written.        
    """

    fin = open(source, 'r')
    num_replaced = count

    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')

    for line in fin:
        out = line
        if pattern in line:
            out = ' '+replace+'\n'
        fout.write(out)

        if out != line:
            num_replaced += 1
        if count and num_replaced > count:
            break
    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E

    fin.close()
    fout.close()

    if not dest:
        shutil.move(name, source)

# ---Module regime consts and variables---
print_prefix='lib.utils>>'

def decode_depth(soil_str):
    strt = int(soil_str[2:5])*0.01
    end = int(soil_str[5:])*0.01
    return strt, end
def fetch_pkgdata(fn_path):
    try:
        data_file = pkg_resources.resource_filename(
            const.PKG_NAME, fn_path)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Config file '{fn_path}' not found in '{const.PKG_NAME}'.")
    return data_file

def valid_path(path):
    if len(path)<5:
        throw_error(f'invalid path: {path}')    
    return path
# ---Classes and Functions---

def throw_error(msg):
    '''
    throw error and exit
    '''
    logging.error(msg)
    exit()

def write_log(msg, lvl=20):
    '''
    write logging log to log file
    level code:
        CRITICAL    50
        ERROR   40
        WARNING 30
        INFO    20
        DEBUG   10
        NOTSET  0
    '''

    logging.log(lvl, msg)

# ---Unit test---
if __name__ == '__main__':
    pass

