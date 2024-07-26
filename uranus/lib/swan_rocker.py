#/usr/bin/env python3
"""Build SWAN workflow"""
import os, shutil, datetime, subprocess
import pandas as pd
import xarray as xr
import numpy as np
import re
import netCDF4 as nc4
import wrf
from . import utils, io, const, cfgparser

# ---Module regime consts and variables---
print_prefix='lib.SWANRocker>>'


# ---Classes and Functions---
class SWANRocker:

    '''
    Construct icbcmaker to generate initial and boundary file for SWAN
    '''
        
    def __init__(self, uranus):
        """ construct maker obj """
        self.uranus=uranus
        self.cfg=uranus.cfg
        self.mach_name=self.cfg['URANUS']['machine_name']
        swancfg=self.cfg['SWAN']
        # central controller
        self.rock_swan=swancfg.getboolean('rock_swan')
        self.run_swan=swancfg.getboolean('run_swan')
        self.restart_run=swancfg.getboolean('restart_run')
        self.swan_root=uranus.cplexe_root
        self.domid=uranus.domain_lv
        self.run_maker=swancfg.getboolean('gen_swan_icbc')
        self.strt_time,self.end_time=uranus.sim_strt_time,uranus.sim_end_time
        self.seglen=float(swancfg['seg_len'])
        
        self.fprefix='swan_bdy'
        self.rewrite_wind=swancfg.getboolean('rewrite_wind')
        self.wrf_dir=utils.parse_fmt_timepath(
            self.strt_time, swancfg['wind_path'])
        self.wind_time_delta=datetime.timedelta(
            minutes=int(swancfg['wind_time_delta']))
        self.domdb_root=uranus.domdb_root
        self.nml_temp=self.cfg['URANUS']['nml_temp']
        self.drv_type=swancfg['drv_type']
        self.drv_dic=const.DRV_DIC[self.drv_type+'_swan']
        self.drv_root=utils.parse_fmt_timepath(
            self.strt_time, swancfg['drv_root'])
        self.proj_root=uranus.proj_root

        self.wavfn_lst, self.file_time_series=io.gen_patternfn_lst(
            self.drv_root, self.drv_dic, self.strt_time, self.end_time, 
            kw='wav',include='both')
        
        self.mach_meta=self.uranus.machine_dic
        # sanity check
        stat=io.check_filelist(self.wavfn_lst)
        if not(stat):
            if self.uranus.fcst_flag:
                for buf_day in range(1,3):
                    utils.write_log(
                        f'{print_prefix}{self.drv_root} missing wav file(s), try {buf_day} days back', 30)
                    buf_time=self.strt_time+datetime.timedelta(days=-buf_day)
                    self.drv_root=utils.parse_fmt_timepath(buf_time, self.cfg['SWAN']['drv_root'])
                    self.wavfn_lst, self.file_time_series=io.gen_patternfn_lst(
                        self.drv_root, self.drv_dic, self.strt_time, self.end_time, kw='wav')
                    stat=io.check_filelist(self.wavfn_lst)
                    if stat:    
                        break
                if not(stat):
                    utils.throw_error(f'{print_prefix}{self.drv_root} missing wav file(s)')
            else:
                utils.throw_error(f'{print_prefix}{self.drv_root} missing wav file(s)')
        
        utils.write_log(f'{print_prefix}SWANMaker Initiation Done.')
    def rock(self):
        """ build initial and boundary files for SWAN """
    
        
        if self.rock_swan:
            if self.run_maker:
                self.clean_workspace()
                # load domain file
                self.load_domain()
                # build segs
                self.build_segs()
                self.print_seg_cmd()
                # parse seg waves from file
                self.parse_seg_waves()
                # generate seg boundary files
                self.gen_seg_files()
                
       
            if self.rewrite_wind and self.uranus.ntasks_atm==0:
                self.interp_wind()

            if self.run_swan:
                self.exeswan()
    
    def exeswan(self):
        '''run swan'''
        mach_meta=self.mach_meta
        bashrc=mach_meta['bashrc']
        mpicmd=mach_meta['mpicmd']
        swan_np=self.uranus.ntasks_wav
        cmd=utils.build_execmd(
            self.mach_name, bashrc, self.swan_root, mpicmd, swan_np, 'coawstM')
        relative_cpl_in=os.path.join(
                '.','Projects',self.nml_temp,'swan_d01.in')
        cmd+=f' {relative_cpl_in}'
        cmd+=f' >& {self.swan_root}/coawstM.log' 
        utils.write_log(print_prefix+'Run swan.exe: '+cmd)
        subprocess.run(cmd, shell=True)   
        
         
    def load_domain(self):
        """ load domain file """
        utils.write_log(print_prefix+'Load domain file...')
        domfn=os.path.join(
                self.domdb_root, self.nml_temp, 'roms_d01_omp.nc')
    
        ds_swan=xr.load_dataset(domfn)
        self.lat2d=ds_swan['lat_rho'].values
        self.lon2d=ds_swan['lon_rho'].values
        self.mask=ds_swan['mask_rho'].values
        
        # resolution in degree
        res_deg=abs(self.lat2d[1,0]-self.lat2d[0,0])
        self.max_seglen=int(self.seglen/res_deg)
        utils.write_log(print_prefix+'Max seg len: %d' % self.max_seglen)
    
    def build_segs(self):
        """ build_segs for SWAN """
        utils.write_log(print_prefix+'build segments...')
        self.segs=[]
        # uid for segs
        self.uid=0
        # 4 boundaries, take 2px width of mask boundary
        self.form_bdy('W', self.mask[:,:2],
                    self.lat2d[:,0], self.lon2d[:,0])
        self.form_bdy('S', self.mask[:2,1:],
                    self.lat2d[0,1:], self.lon2d[0,1:])
        self.form_bdy('E', self.mask[1:,-2:],
                    self.lat2d[1:,-1], self.lon2d[1:,-1])
        self.form_bdy('N', self.mask[-2:,1:-2],
                    self.lat2d[-1,1:-2], self.lon2d[-1,1:-2])
        
        for seg in self.segs:
            seg['file']=self.fprefix+'.%s.%03d.txt' % (seg['orient'], seg['id']) 



    def prepare_cplrock(self):
        domfn=os.path.join(
                self.uranus.domdb_root, self.uranus.nml_temp, 'swan*')
        io.copy_files(domfn, self.proj_root)
        
        nml_src=os.path.join(
                self.uranus.cfgdb_root, self.uranus.nml_temp,f'swan_{self.domid}.in')
        nml_dest=os.path.join(self.proj_root, f'swan_{self.domid}.in')
        shutil.copy(nml_src, nml_dest)
        
        utils.sedline(
            'ssyyyymmdd.hh',self.strt_time.strftime('%Y%m%d.%H'),
            nml_dest,whole_line=False)
        utils.sedline(
            'eeyyyymmdd.hh',self.end_time.strftime('%Y%m%d.%H'),
            nml_dest,whole_line=False)
        if self.restart_run:
            utils.sedline(
                '&INITIAL','INITIAL',nml_dest,whole_line=False)
    def clean_workspace(self):
        io.del_files(self.swan_root, const.SWAN_CLEAN_LIST)
            
            
    def print_seg_cmd(self):
        """ print seg cmd for swan.in 
        """
        utils.write_log(print_prefix+'print seg cmd for swan.in...')
        cmd_line=''
        for seg in self.segs:
            cmd_line='%sBOUNDSPEC SEGMENT XY %8.4f %8.4f %8.4f %8.4f VARIABLE FILE 0 \'%s\'\n' % (
                cmd_line, seg['lon0'], seg['lat0'], seg['lon1'], seg['lat1'], seg['file'])
        
        with open(self.proj_root+'/swan_d01.in', 'r') as sources:
            lines = sources.readlines()

        with open(self.proj_root+'/swan_d01.in', 'w') as sources:
            for line in lines:
                # regexp pipeline
                line=re.sub('@BOUNDSPEC', cmd_line, line)
                sources.write(line)
    
    def form_bdy(self, bdy_type, maskline, latline, lonline):
        """ form boundary accourding to maskline """
        find_flag=False
        uid=self.uid
        if maskline.shape[0]==2:
            maskline=maskline.T
        for i in range(maskline.shape[0]):
            if (maskline[i,0] == 1) and (maskline[i,1]==1): # ocean point
                if not(find_flag):
                    find_flag=True
                    seg_dict={'id':uid, 'orient':bdy_type, 
                    'lat0':latline[i], 'lon0':lonline[i]}
                    uid=uid+1
                    seg_len=1
                else:
                    seg_len=seg_len+1
                    if seg_len==self.max_seglen:
                        seg_dict=close_seg(seg_dict, 
                            latline[i], lonline[i], seg_len)
                        self.segs.append(seg_dict)
                        seg_dict={}
                        find_flag=False
            else: # find land point on the boundary (width=2px)
                # already in seg
                if find_flag:
                    if seg_len>int(0.25*self.max_seglen):
                        seg_dict=close_seg(seg_dict, 
                            latline[i-1], lonline[i-1], seg_len)
                        self.segs.append(seg_dict)
                    seg_dict={}
                    find_flag=False
            # last position
            if i==maskline.shape[0]-1:
                if find_flag:
                    seg_dict=close_seg(seg_dict,
                        latline[i], lonline[i], seg_len)
                    self.segs.append(seg_dict)
                    seg_dict={}

    
    
    def parse_seg_waves(self):
        """ parse seg cmd for swan.in 
        """
        utils.write_log(print_prefix+'parse seg waves from bdy files...')
               # read the first file
        
        ds_grib=[] 
        for fn in self.wavfn_lst: 
           ds_grib.append(xr.load_dataset(fn, engine='cfgrib', backend_kwargs={'errors': 'ignore'}))
        
        comb_ds=xr.concat(ds_grib, 'time')
        
        src_lon=comb_ds['longitude'].values
        src_lat=comb_ds['latitude'].values
        
        utils.write_log(print_prefix+'swan lat min: %f, max: %f; lon min: %f, max: %f' % (
            np.min(self.lat2d), np.max(self.lat2d), np.min(self.lon2d), np.max(self.lon2d)))
        utils.write_log(print_prefix+'bdy lat min: %f, max: %f; lon min: %f, max: %f' % (
            np.min(src_lat), np.max(src_lat), np.min(src_lon), np.max(src_lon)))
 
        if not(is_domain_within_bdyfile(self.lat2d, self.lon2d, src_lat, src_lon)):
            utils.throw_error(print_prefix+'SWAN domain is out of wave bdy file! Cannot continue!')
        src_lat2d,src_lon2d=np.meshgrid(src_lat, src_lon)
        src_mask=np.isnan(comb_ds['swh'].isel(time=0).values)
        ts=comb_ds['valid_time'].values
        len_ts=len(ts)
        for seg in self.segs:
            i,j=find_ij_mask(src_lat2d,src_lon2d,src_mask,seg['latm'],seg['lonm'])
            seg['sigh']=comb_ds['swh'].values[:,j,i]
            seg['sigh']=np.where(seg['sigh']<0.1, 0.1, seg['sigh'])
            if self.drv_type=='gfs':
                seg['period']=comb_ds['perpw'].values[:,j,i]
                seg['dir']=comb_ds['dirpw'].values[:,j,i]
                seg['spread']=np.repeat(20.0, len_ts)
            else:
                seg['period']=comb_ds['mwp'].values[:,j,i]
                seg['dir']=comb_ds['mwd'].values[:,j,i]
                seg['spread']=cal_dir_spread(comb_ds['wdw'].values[:,j,i])
            data=np.array([seg['sigh'], seg['period'], seg['dir'], seg['spread']])
            seg['df'] = pd.DataFrame(data.T, index=ts, columns=['sigh', 'period', 'dir', 'spread'])
           
    def gen_seg_files(self):
        """ generate seg files """ 
        utils.write_log(print_prefix+'generate seg files...')
        
        for seg in self.segs:
            seg_file=open(self.swan_root+'/'+seg['file'],'w')
            #seg_file=open('test.csv','w')
            seg_file.write('TPAR\n')
            for tp in seg['df'].itertuples():
                seg_file.write('%s %8.2f %8.2f %8.2f %8.2f\n' 
                % (tp[0].strftime('%Y%m%d.%H%M'), tp[1], tp[2], tp[3], tp[4]))
            seg_file.close()
    
    
    def interp_wind(self):
        """ interpolate wind from WRFOUT to SWAN """
        utils.write_log(print_prefix+'Interpolate wind for SWAN...')
        cfg=self.cfg 
        # WRF Parameters
        wrf_dir=self.wrf_dir
        if not os.path.exists(wrf_dir):
            utils.throw_error('WRF output directory does not exist!')
        
        wrf_domain=cfg['SWAN']['wrf_match']
        wind_prefix=cfg['SWAN']['wind_prefix']
        domdb_path=os.path.join(self.uranus.domdb_root, self.uranus.nml_temp) 
        # SWAN domain file
        dom_file=os.path.join(domdb_path,f'roms_{self.domid}_omp.nc')
        if not(os.path.exists(dom_file)):
            dom_file=os.path.join(domdb_path,f'swan_{self.domid}.nc')
        ds_swan=xr.load_dataset(dom_file)

        lat_swan=ds_swan['lat_rho']
        lon_swan=ds_swan['lon_rho']
        
        force_fn=wind_prefix+'_'+self.domid+'.dat'
        force_fn=self.proj_root+'/'+force_fn
        
        # IF force file exists
        if os.path.exists(force_fn):
            utils.write_log(print_prefix+'Delete existing wind file...%s' % force_fn, 30)
            os.remove(force_fn)
        
        curr_time=self.strt_time

        # iter timeframes 
        while curr_time < self.end_time:
            wrf_file=io.get_wrf_fn(curr_time, wrf_domain)
            wrf_file=os.path.join(wrf_dir,wrf_file)
            utils.write_log('Read '+wrf_file)
            
            utils.write_log(print_prefix+'Read WRFOUT surface wind...')
            wrf_hdl=nc4.Dataset(wrf_file)
            
            wrf_u10 = wrf.getvar(
                    wrf_hdl, 'U10', 
                    timeidx=wrf.ALL_TIMES, method="cat")
            wrf_v10 = wrf.getvar(
                    wrf_hdl, 'V10',
                    timeidx=wrf.ALL_TIMES, method="cat")
            wrf_time=wrf.extract_times(
                    wrf_hdl,timeidx=wrf.ALL_TIMES, do_xtime=False)
            wrf_hdl.close()
            
            wrf_time=[pd.to_datetime(itm) for itm in wrf_time]
            
            utils.write_log(print_prefix+'Query Wind Timestamp:'+str(curr_time))
            u10_tmp=wrf_u10
            v10_tmp=wrf_v10

            swan_u = interp_wrf2swan(u10_tmp, lat_swan, lon_swan)
            swan_v = interp_wrf2swan(v10_tmp, lat_swan, lon_swan)
            
            if 'swan_uv' in locals():# already defined
                swan_uv=np.concatenate((swan_uv, swan_u, swan_v), axis=0)
            else:
                swan_uv=np.concatenate((swan_u, swan_v), axis=0)
            
            curr_time=curr_time+self.wind_time_delta
            
        # add one more time stamp to close up the file
        swan_uv=np.concatenate((swan_uv, swan_u, swan_v), axis=0)
                
        utils.write_log('Output...')
        with open(force_fn, 'a') as f:
            np.savetxt(f, swan_uv, fmt='%7.2f', delimiter=' ')
            f.write('\n')
        del swan_uv
def is_domain_within_bdyfile(lat_swan, lon_swan, lat_bdy, lon_bdy):
    '''
    test if swan domain is within the boundary file domain
    '''
   
    if not(lat_swan.min()>=lat_bdy.min() and lat_swan.max()<=lat_bdy.max()):
        write_log('lat_swan:%7.4f-%7.4f is not within lat_bdy:%7.4f-%7.4f'%(
            lat_swan.min(), lat_swan.max(), lat_bdy.min(), lat_bdy.max()))
        return False
    if not(lon_swan.min()>=lon_bdy.min() and lon_swan.max()<=lon_bdy.max()):
        return False
    return True

def find_ij_mask(lat2d, lon2d, mask, latm, lonm):
    """ find i,j in src_lat, src_lon, src_mask 
    """
    dislat=lat2d-latm
    dislon=lon2d-lonm
    dis=np.sqrt(dislat**2+dislon**2)
    ids=np.argsort(dis, axis=None)
    ij_tp=np.unravel_index(ids,lat2d.shape)
    for i,j in zip(ij_tp[0], ij_tp[1]):
        # over ocean
        if mask[j,i]==False:
            #utils.write_log(print_prefix+'search latm=%9.2f,lonm=%9.2f...' % (latm, lonm))
            #utils.write_log(print_prefix+'find i=%3d,j=%3d...'%(i,j))
            #utils.write_log(print_prefix+'find latx=%9.2f,lonx=%9.2f, with dis=%9.2f' 
            #    % (lat2d[i,j], lon2d[i,j], dis[i,j]*111))
            return i,j

def cal_dir_spread(sigma):
    '''
    return parameterized spread angle 
    '''
    return  np.arcsin(sigma/np.sqrt(2))*180.0/np.pi

           
def close_seg(seg_dict, lat, lon, seg_len):
    """close a segment"""
    seg_dict['lat1']=lat
    seg_dict['lon1']=lon
    seg_dict['latm']=(seg_dict['lat0']+seg_dict['lat1'])/2
    seg_dict['lonm']=(seg_dict['lon0']+seg_dict['lon1'])/2
    seg_dict['len']=seg_len
    return seg_dict



def interp_wrf2swan(wrf_var, swan_lat, swan_lon):
    """ 
    Linearly interpolate var from WRF grid onto SWAN grid 
    """
    
    from scipy import interpolate

    x_org=wrf_var.XLAT.values.flatten()
    y_org=wrf_var.XLONG.values.flatten()
    
    interp=interpolate.LinearNDInterpolator(list(zip(x_org, y_org)), wrf_var.values.flatten())
    #interp=interpolate.NearestNDInterpolator(list(zip(x_org, y_org)), wrf_var.values.flatten())
    template = interp(swan_lat.values, swan_lon.values)
    return template  