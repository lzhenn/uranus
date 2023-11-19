import numpy as np
from . import utils, const

def sigma2depth(zeta,h,ds_smp):
        
    # S-coordinate parameter, critical depth
    sc = ds_smp.sc_r.values
    Cs = ds_smp.Cs_r.values
    hc = ds_smp.hc.values
    vtrans=ds_smp.Vtransform.values

    nz=len(Cs)
    nt,nx,ny=zeta.shape
    
    if vtrans == 1:
        Zo_rho = hc * (sc - Cs) + Cs * h 
        z_rho = Zo_rho + zeta * (1 + Zo_rho / h)
    elif vtrans == 2:
        h_3d=np.broadcast_to(h, (nz, nx, ny))
        zeta_3d=np.broadcast_to(zeta, (nz, nx, ny))
        Zo_rho =np.zeros((nz,nx,ny)) 
        for i in range(len(Cs)):
            Zo_rho[i,:,:]=(hc * sc[i] + Cs[i] * h_3d[i,:,:]) / (hc + h_3d[i,:,:])
            z_rho = zeta_3d + Zo_rho * (zeta_3d + h_3d)
    else:
        utils.throw_error('vtrans must be 1 or 2')
        
    return z_rho

def assign_depth_idx(z_rho, mask):
    '''
    Assign depth index to each grid point with super fast method
    '''
    
    nz,nx,ny=z_rho.shape
    NZ_DP=len(const.DENSE_DP)
    mask3d=np.broadcast_to(mask, (nz, nx, ny))
    z_rho=(-1.0*z_rho).astype(int)
    dp_idx=z_rho.copy()
    
    # ----super fast assign depth index----
    # First remove land points 
    dp_idx=np.where(mask3d,dp_idx,0)
    # Then assign depth index, is this clear why we keep a constant DENSE_DP?
    dp_idx=np.where(z_rho>50,45+z_rho//10,dp_idx)
    dp_idx=np.where(z_rho>200,63+z_rho//100,dp_idx)
    dp_idx=np.where(z_rho>1000,71+z_rho//500,dp_idx)
    # upper bound
    dp_idx=np.where(dp_idx>=NZ_DP,NZ_DP-1,dp_idx)
    return dp_idx


def hybrid2pressure(da,ap,b,ps,PLVS):
    '''
    Convert hybrid level to pressure level, fast nearest interpolation
    '''
    #print(da.values.shape)
    nz, nlat, nlon=da.values.shape
    da_new=da.copy()[0:len(PLVS),:,:]
    try:
        da_new=da_new.rename({'lev':'plev'})
    except ValueError:
        pass
    da_new=da_new.assign_coords({'plev':PLVS})
    pa3d=np.broadcast_to(ps.values, (nz, nlat, nlon))
    pa3d=np.swapaxes(pa3d,0,2) 
    # get hybrid level
    pa3d=pa3d*b+ap
    azflag=False
    if pa3d[0,0,0]<pa3d[0,0,nz-1]:
        azflag=True
    for idz, plv in enumerate(PLVS):
        idx2d=np.sum(np.where(pa3d<plv,0,1),axis=2)-1
        idx2d=np.where(idx2d<0,0,idx2d)
        if azflag: # small to big
            idx2d=nz-idx2d-1
        idx3d=np.broadcast_to(idx2d.T, (nz, nlat, nlon))
        temp=np.take_along_axis(da.values, idx3d, axis=0)
        da_new.values[idz,:,:]=temp[0,:,:]
    return da_new


def find_indices(strt_dp, end_dp, SOI_LVS):
    # Check if strt_dp is 0
    if strt_dp == 0:
        i = 0
        res_s=0
    else:
        i = next(x[0] for x in enumerate(SOI_LVS) if x[1] >= strt_dp)
        res_s=strt_dp-SOI_LVS[i-1]

    # Check if end_dp is greater than max value in SOI_LVS
    if end_dp > max(SOI_LVS):
        j = len(SOI_LVS) - 1
        res_e=0
    else:
        j = next(x[0] for x in enumerate(SOI_LVS) if x[1] >= end_dp)
        res_e=SOI_LVS[j]-end_dp
    return i, j, res_s, res_e
