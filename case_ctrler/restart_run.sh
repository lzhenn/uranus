#/bin/sh 
ID=$1
PID=$2
HOST=hqlx$ID
# Set the path to your wrfrst files
wrf_path=/home/lzhenn/array${ID}/WRF-4.3_P${PID}/run

# Get the latest wrfrst file
latest_file=$(ls -t $wrf_path/wrfrst_d04_* | head -n1)
latest_file=`basename $latest_file`
# Extract the latest data from the file
latest_date=${latest_file:11:10}

echo "The latest wrfrst data in $wrf_path is: $latest_date"

# Extract the year, month, and day from the new date
new_year=$(date -d $latest_date +%Y)
new_month=$(date -d $latest_date +%m)
new_day=$(date -d $latest_date +%d)

namelist_path=$wrf_path/namelist.input
# Use sed to update the start date entries in the namelist.input file
sed -i "s/start_year.*/start_year = $new_year,$new_year,$new_year,$new_year/" $namelist_path
sed -i "s/start_month.*/start_month = $new_month,$new_month,$new_month,$new_month/" $namelist_path
sed -i "s/start_day.*/start_day = $new_day,$new_day,$new_day,$new_day/" $namelist_path
sed -i "s/restart.*.false./restart = .true./" $namelist_path


ssh $HOST "source /home/lzhenn/.bashrc_intel20_amd;cd ${wrf_path}; mpirun -np 64 ./wrf.exe"