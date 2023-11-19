import uranus
#uranus.copy_cfg('config.cmip6_wrf.p2.ini') # for test 
agent=uranus.Uranus(cfgfn='config.cmip6_wrf.111p1.ini')
agent.waterfall()
