import uranus
#uranus.copy_cfg('config.cpsv3.ini') # for test 
agent=uranus.Uranus(cfgfn='config.cpsv3.ini')
agent.waterfall()
