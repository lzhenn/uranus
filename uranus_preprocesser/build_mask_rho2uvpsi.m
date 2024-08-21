
%addpath(genpath('/home/metctm1/array/project/1911-COAWST/script/mfiles'));


roms_grid='/home/lzhenn/array74/Njord_Calypso/COAWST_Poseidon/Projects/Poseidon/roms_d03.org.nc'; % Call generic grid creation.

netcdf_load(roms_grid)
% Compute mask on rho, u, v, and psi points

water = double(mask_rho);
u_mask = water(1:end-1,:) & water(2:end,:);
v_mask= water(:,1:end-1) & water(:,2:end);
psi_mask= water(1:end-1,1:end-1) & water(1:end-1,2:end) & water(2:end,1:end-1) & water(2:end,2:end);
ncwrite(roms_grid,'mask_rho',mask_rho);
ncwrite(roms_grid,'mask_u',double(u_mask));
ncwrite(roms_grid,'mask_v',double(v_mask));
ncwrite(roms_grid,'mask_psi',double(psi_mask));

