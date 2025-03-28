import uranus
agent=uranus.Uranus(cfgfn='config.SEA.hq206.ini')
#agent=uranus.Uranus(cfgfn='config.SEA.liquor.preprocess.ini')
#agent=uranus.Uranus(cfgfn='config.SEA.hq84.ini')
#agent=uranus.Uranus(cfgfn='config.SEA.hq182.ini')
#agent=uranus.Uranus(cfgfn='config.SEA.monpack.ini')
agent.waterfall()
