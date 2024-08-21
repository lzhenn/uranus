#/usr/bin/env python3
'''
Date: June 16, 2021
Initial configuration for fundamental parameters
Zhenning LI
'''

import os
import lib.cfgparser

sys_cfg_fn='./glb_cfg/config_sys.ini'
cfg_sys=lib.cfgparser.read_cfg(sys_cfg_fn)

# write root 
cfg_sys['SYS']['root']=os.getcwd()
lib.cfgparser.write_cfg(cfg_sys, sys_cfg_fn)


