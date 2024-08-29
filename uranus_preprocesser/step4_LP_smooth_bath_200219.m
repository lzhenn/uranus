%addpath(genpath('/home/metctm1/array/project/1911-COAWST/script/mfiles'));
%raw_fn='/home/lzhenn/Njord_dev/Projects/Njord_t123/roms_d03.nc';
%out_fn='/home/lzhenn/Njord_dev/Projects/Njord_t123/roms_d03.lp0d1.nc';
%raw_fn='/home/lzhenn/array74/workspace/uranus/uranus/domaindb/njord_noluzon/roms_d01.nc';
%out_fn='/home/lzhenn/array74/workspace/uranus/uranus/domaindb/njord_noluzon/roms_d01_omp.nc';
raw_fn='/home/lzhenn/array74/workspace/uranus/uranus/domaindb/les_aqua/roms_d01.nc';
out_fn='/home/lzhenn/array74/workspace/uranus/uranus/domaindb/les_aqua/roms_d01_omp.nc';

[status,cmdout] = system(['cp ' raw_fn ' ' out_fn]);
status

netcdf_load(raw_fn);

rx0max=0.1;
disp(['Target for rx0=' num2str(rx0max)]);

disp('Using LP method with heuristic');

%h0=GRID_SmoothPositive_rx0(mask_rho, h, rx0max);
h0=GRID_LinProgHeuristic_rx0(mask_rho, h, rx0max);

ncwrite(out_fn,'h',h0);
figure
pcolorjw(lon_rho,lat_rho,h0-h)
hold on
caxis([-10 10]); colorbar
title('ROMS bathy diff')
xlabel('longitude'); ylabel('latitiude')
