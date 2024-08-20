import uranus
import argparse
import configparser

config=configparser.ConfigParser()
config.read('config.yt2511.ini')
agent=uranus.Uranus(cfg=config)

agent.waterfall()
