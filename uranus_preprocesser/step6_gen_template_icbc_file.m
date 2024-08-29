% roms_master_climatology_coawst_mw
%
% This routine :
%  - creates climatology, boundary, and initial condition files for ROMS: 
%    coawst_clm.nc ; coawst_bdy.nc ; coawst_ini.nc 
%    on a user-defined grid for a user-defined date.
%
% This is currently set up to use opendap calls to acquire data
% from HYCOM + NCODA Global 1/12 Degree Analysis and interp to roms grid.
%  
% based on efforts by:
% written by Mingkui Li, May 2008
% Modified by Brandy Armstrong March 2009
% jcwarner April 20, 2009
% Ilgar Safak modified on June 27, 2012 such that now:
% - HYCOM url is a user-definition
% - "hc" is called from the structure "gn".(still needs to be tested with wet/dry).
% - updatinit_coawst_mw.m modified to get desired time (T1) as a variable;
%    ocean_time=T1-datenum(1858,11,17,0,0,0)
% Updates from Christie Hegermiller, Feb 2019
%

%%%%%%%%%%%%%%%%%%%%%   START OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
domain_str='d01';


% (1) Enter start date (T1) and number of days to get climatology data 
%T1 = datetime(2021,06,15,12,0,0);
T1=datetime(2017,05,23,00,0,0);
%number of days and frequency to create climatology files for
%numdays = 3;
numdays=1;
dayFrequency = 1;

% (2) Enter URL of the HYCOM catalog for the requested time, T1
%     see http://tds.hycom.org/thredds/catalog.html
%url = '/users/b145872/project-dir/data/hycom/goni/';
%url = '/home/lzhenn/array74/Njord_Calypso//drv_field//hycom_subset/'; 
ocn_ra= '/home/lzhenn/drv_field//hycom_subset/2017052300/';

% (3) Enter working directory (wdr)
%wdr = ['/home/metctm1/array/app/COAWST/COAWST_operational/Projects/GBA_operational/ow_icbc/', domain_str];
wdr = '/home/lzhenn/array74/Njord_Calypso/domaindb/les_aqua/';
%roms_swan_grid_dir='/users/b145872/project-dir/app/COAWST-FULL/Projects/GONI/grid/';
%roms_swan_grid_dir='/home/metctm1/array/app/COAWST/COAWST_operational/Projects/GBA_operational/roms_swan_grid/';
roms_swan_grid_dir = '/home/lzhenn/array74/Njord_Calypso/domaindb/les_aqua/';

% (4) Enter path and name of the ROMS grid
%modelgrid = [roms_swan_grid_dir,'roms_',domain_str,'_lp0d1.nc'];
modelgrid = [roms_swan_grid_dir,'roms_d01_omp.nc'];
disp(['ROMS grid: ', modelgrid])
% (5) Enter grid vertical coordinate parameters --These need to be consistent with the ROMS setup. 
theta_s     =  6.0;
theta_b     =  0.4;
Tcline      = 50.0;
N           =  30;
Vtransform  =  2;
Vstretching =  4;

%%%%%%%%%%%%%%%%%%%%%   END OF USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%

% deal with swan grid first
disp('generate swan grid first...')
eval(['cd ',roms_swan_grid_dir])
%roms2swan(modelgrid)

%roms icbc
eval(['cd ',wdr])

tic

% Call to get HYCOM indices for the defined ROMS grid
disp('getting roms grid, hycom grid, and overlapping indices')
[gn, clm]=get_ijrg_exp930([ocn_ra,'hycom.grid.nc'], modelgrid, theta_s, theta_b, Tcline, N, Vtransform, Vstretching);

% Call to create the climatology (clm) file
disp('going to create clm file')
fn=updatclim_coawst_mw_local_exp930(T1, gn, clm, 'roms_d01_clmsmp.nc', wdr, ocn_ra)

% Call to create the boundary (bdy) file
disp('going to create bndry file')
updatbdry_coawst_mw(fn, gn, 'roms_d01_bdysmp.nc', wdr)

% Call to create the initial (ini) file
disp('going to create init file')
updatinit_coawst_mw(fn, gn,'roms_d01_inismp.nc', wdr, datenum(T1))

toc


%% Call to create the long climatology (clm) file
if numdays>1
    disp('going to create more days of clm and bnd files')
    if (ispc)
      eval(['!copy coawst_clm.nc coawst_clm_',datestr(T1,'yyyymmdd'),'.nc'])
      eval(['!copy coawst_bdy.nc coawst_bdy_',datestr(T1,'yyyymmdd'),'.nc'])
    else
      eval(['!cp coawst_clm.nc coawst_clm_',datestr(T1,'yyyymmdd'),'.nc'])
      eval(['!cp coawst_bdy.nc coawst_bdy_',datestr(T1,'yyyymmdd'),'.nc'])
    end
    for it=dayFrequency:dayFrequency:numdays-1      %1st day already created, NEED to set number of days at top!
        fname=['coawst_clm_',datestr(T1+it,'yyyymmdd'),'.nc']
        fn=updatclim_coawst_mw_local_exp930(T1+it,gn,clm,fname,wdr,ocn_ra)
        fname=['coawst_bdy_',datestr(T1+it,'yyyymmdd'),'.nc'];
        updatbdry_coawst_mw(fn,gn,fname,wdr)
    end
end

toc
exit

