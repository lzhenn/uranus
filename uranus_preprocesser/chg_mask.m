fn='/home/lzhenn/array74/Njord_Calypso/domaindb/poseidon_LTtmr_L12/roms_d03.nc'

% Open the ROMS netCDF file for writing
ncid = netcdf.open(fn, 'WRITE');

% Read the "mask_rho" variable
mask_rho = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'mask_rho'));

% Compute the "u_mask" variable
u_mask = mask_rho(1:end-1,:) & mask_rho(2:end,:);
% Compute the "v_mask" variable
v_mask = mask_rho(:,1:end-1) & mask_rho(:,2:end);
% Compute the "psi_mask" variable
psi_mask = mask_rho(1:end-1,1:end-1) & mask_rho(1:end-1,2:end) & mask_rho(2:end,1:end-1) & mask_rho(2:end,2:end);


ncwrite(fn,'mask_u',double(u_mask));
ncwrite(fn,'mask_v',double(v_mask));
ncwrite(fn,'mask_psi',double(psi_mask));


% Close the netCDF file
netcdf.close(ncid);
