import uranus
uranus.copy_cfg('config.poseidon.ini') # for test 
agent=uranus.Uranus(cfgfn='config.poseidon.ini')
agent.waterfall()
