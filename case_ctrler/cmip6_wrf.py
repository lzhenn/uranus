import uranus
import argparse
import configparser

config=configparser.ConfigParser()
config.read('config.cmip6_wrf.ini')
print(config)
# Create an ArgumentParser object
parser = argparse.ArgumentParser()
# Add the command-line arguments
parser.add_argument("-m", "--mach", type=str,
                    help="machine ID")
parser.add_argument("-p", "--parallel", type=str,
                    help="parallel ID")
parser.add_argument("-y", "--year", type=str,
                    help="sim year")
parser.add_argument("-d", "--day", type=str,
                    help="sim days")
parser.add_argument("-s", "--scenario", type=str,
                    help="sim scenario")

# Parse the command-line arguments
args = parser.parse_args()
# Print the values of the options
print(args)


config['URANUS']['machine_name']=f'hqlx{args.mach}'
config['URANUS']['model_init_ts']=f'{args.year}010100'
config['URANUS']['model_run_days']=f'{args.day}'
if args.parallel=='1':
    config['WRF']['wps_root']=f'/home/lzhenn/array{args.mach}/WPS-4.3/'
    config['WRF']['wrf_root']=f'/home/lzhenn/array{args.mach}/WRF-4.3/run'
else:
    config['WRF']['wps_root']=f'/home/lzhenn/array{args.mach}/WPS-4.3_P{args.parallel}/'
    config['WRF']['wrf_root']=f'/home/lzhenn/array{args.mach}/WRF-4.3_P{args.parallel}/run'
config['bcmm']['scenario_name']=args.scenario
# Print the contents of the config file
for section in config.sections():
    print(f"[{section}]")
    for option in config.options(section):
        print(f"{option} = {config.get(section, option)}")

agent=uranus.Uranus(cfg=config)

agent.waterfall()
