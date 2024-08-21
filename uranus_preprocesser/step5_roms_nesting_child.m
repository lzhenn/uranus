%addpath(genpath('/home/metctm1/array/project/1911-COAWST/script/mfiles'));
Jstr=136; Jend=244; Istr=260; Iend=386;
ref_ratio=5;
out_root='/home/lzhenn/array74/Njord_Calypso/COAWST_Njord_dev/Projects/Njord_t123/'
roms_child_grid=[out_root,'roms_d03.nc'];
roms_out_grid='/home/lzhenn/Njord_dev/Projects/Njord_t1t2/roms_swan_grid/roms_d02_lp0d1.nc'
contact_name=[out_root,'roms_d0203_contact.nc'];

F=coarse2fine(roms_out_grid,roms_child_grid, ...
              ref_ratio,Istr,Iend,Jstr,Jend);
Gnames={roms_out_grid,roms_child_grid};
[S,G]=contact(Gnames,contact_name);
