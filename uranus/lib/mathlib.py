

def hybrid2pressure(da,ap,b,ps):
    '''
    Convert hybrid level to pressure level
    '''
    #print(da.values.shape)
    nz, nlat, nlon=da.values.shape
    da_new=da.copy()[0:len(PLVS),:,:]
    da_new=da_new.rename({'lev':'plev'})
    da_new=da_new.assign_coords({'plev':new_lv})
    
    pa3d=np.broadcast_to(ps.values, (nz, nlat, nlon))
    pa3d=np.swapaxes(pa3d,0,2) 
    # get hybrid level
    pa3d=pa3d*b+ap
    for idz, plv in enumerate(PLVS):
        idx2d=np.sum(np.where(pa3d<plv,0,1),axis=2)-1
        idx2d=np.where(idx2d<0,0,idx2d)
        idx3d=np.broadcast_to(idx2d.T, (nz, nlat, nlon))
        temp=np.take_along_axis(da.values, idx3d, axis=0)
        da_new.values[idz,:,:]=temp[0,:,:]
    return da_new