#/usr/bin/env python3
"""
    Commonly used utilities

    Function    
    ---------------
    throw_error(msg):
        throw error and exit
    
    write_log(msg, lvl=20):
        write logging log to log file
    
    parse_tswildcard(tgt_time, wildcard):
        parse string with timestamp wildcard 
        to datetime object

"""
# ---imports---
import logging
import pkg_resources
from . import const


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

def parse_tswildcard(tgt_time, wildcard):
    '''
    parse string with timestamp wildcard to datetime object
    '''
    seg_str=wildcard.split('@')
    parsed_str=''
    for seg in seg_str:
        if seg.startswith('%'):
            parsed_str+=tgt_time.strftime(seg)
        else:
            parsed_str+=seg
    return parsed_str


# ---Unit test---
if __name__ == '__main__':
    pass

