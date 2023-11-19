# conda create -n uranus86 python=3.10
# conda activate uranus86
# conda install numpy scipy xarray cftime h5netcdf pandas

pip uninstall uranus 
rm -rf build dist *egg-info
python setup.py sdist
pip install .
#cd test_case
#python test_cpsv3_th2.py 
#cd ..
