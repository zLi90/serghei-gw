"""
    Plot SERGHEI-SW output for the heterogeneous slab test problem
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

fdir = 'out/'
ind = [3]
dim = [100, 1, 50]
N = int(dim[0]*dim[1]*dim[2])
jy = -1
zmax = 15.0
zmin = 0.0
dz = 0.1
slope = 0.1


porosity = 0.1

h_out = []
wc_out = []
fname = fdir + 'output_subsurface.nc'
fid = netcdf.NetCDFFile(fname,'r')
print(fid.variables)
h_out = fid.variables['hd']
wc_out = fid.variables['wc']
fid.close()

print(np.shape(h_out))

#   Create z coordinates
xx = np.linspace(0.0, 10.0, dim[0])
zz = np.linspace(zmin, zmax, int((zmax-zmin)/dz))
satu = np.nan*np.ones((len(zz), len(xx)), dtype=float)


"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""

plt.figure(1, figsize=[12, 5])

#   Plot saturation
ifig = 1
for ii in range(len(ind)):
    plt.subplot(1,len(ind),ifig)
    slice = wc_out[ind[ii],:,jy,:]/porosity
    for jj in range(dim[0]):
        for kk in range(dim[2]):
            kk2 = int(kk + slope*(dim[0]-jj)*10.0)
            satu[kk2,jj] = slice[kk,jj]
    plt.imshow(satu, vmin=0.2, vmax=1.0, cmap='jet_r')

    plt.colorbar()
    plt.title('Saturation')
    ifig += 1

#   Plot pressure head
plt.figure(2, figsize=[12, 5])
ifig = 1
for ii in range(len(ind)):
    plt.subplot(1,len(ind),ifig)
    slice = h_out[ind[ii],:,jy,:]
    for jj in range(dim[0]):
        for kk in range(dim[2]):
            kk2 = int(kk + slope*(dim[0]-jj)*10.0)
            satu[kk2,jj] = slice[kk,jj]
    plt.imshow(satu, vmin=-2.0, vmax=2.0, cmap='jet')

    plt.colorbar()
    plt.title('H')
    ifig += 1

plt.show()
