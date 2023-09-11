import numpy as np

def hybrid2pressure(da,ap,b,ps,PLVS):
    '''
    Convert hybrid level to pressure level, fast nearest interpolation
    '''
    #print(da.values.shape)
    nz, nlat, nlon=da.values.shape
    da_new=da.copy()[0:len(PLVS),:,:]
    da_new=da_new.rename({'lev':'plev'})
    da_new=da_new.assign_coords({'plev':PLVS})
    
    pa3d=np.broadcast_to(ps.values, (nz, nlat, nlon))
    pa3d=np.swapaxes(pa3d,0,2) 
    # get hybrid level
    pa3d=pa3d*b+ap
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