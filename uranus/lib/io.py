#!/usr/bin/env python3
"""specific module for IO"""
# ---imports---
import os
import pandas as pd
from . import utils

# ---Module regime consts and variables---
print_prefix='lib.io>>'


# ---Classes and Functions---
def del_files(tgt_path, fnpatterns):
    file_list = os.listdir(tgt_path)
    utils.write_log(f'{print_prefix}Clean workspace for {tgt_path}...')
    for filename in file_list:
        if filename.startswith(tuple(fnpatterns)):
            os.remove(os.path.join(tgt_path, filename))
def gen_patternfn_lst(drv_root, drv_dic, inittime, endtime, kw='atm'):
    fn_lst=[]
    tfs=pd.date_range(start=inittime, end=endtime, freq=drv_dic['atm_file_nfrq'])
    pattern=drv_dic[f'{kw}_naming_pattern']
    if '$I' in pattern:
        init_fmt=pattern.split('$I')[1]
        init_str = inittime.strftime(init_fmt)
        pattern=pattern.replace(f'$I{init_fmt}$I', init_str)
    frm_fmt=pattern.split('$F')[1]
    
    for tf in tfs[:-1]: # do not include the last time frame
        frm_str = tf.strftime(frm_fmt)
        # Replace the placeholders with the formatted strings
        fn=pattern.replace(f'$F{frm_fmt}$F', frm_str)
        fn_full=os.path.join(drv_root, fn)
        if not(os.path.exists(fn_full)):
            utils.throw_error(f'file not exist: {fn_full}')
        fn_lst.append(fn_full)

    return fn_lst
# ---Unit test---
if __name__ == '__main__':
    pass