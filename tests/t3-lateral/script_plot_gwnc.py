"""
    Verify serghei subsurface solver with Tracy's analytical solution

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

fdir = 'out/'
ind = ['0','3','7','10']
flabel = ['0','1','2','3']
jy = -1

dim = [40, 1, 76]
N = int(dim[0]*dim[1]*dim[2])
zVec = np.linspace(-15, 4.0, N)

fname = fdir + 'output_subsurface.nc'
fid = netcdf.NetCDFFile(fname,'r')
h_out = fid.variables['hd']
wc_out = fid.variables['wc']
fid.close()
print(' Shape of the output fields: ', np.shape(h_out))

"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""

plt.figure(1, figsize=[8, 5])

ifig = 1
for ii in range(len(ind)):
    plt.subplot(2,len(ind),ifig)
    plt.imshow(np.transpose(np.transpose(h_out[ii,:,jy,:])), vmin=-10.0, vmax=10.0, cmap='jet')
    plt.colorbar()
    plt.title('H')
    ifig += 1
for ii in range(len(ind)):
    plt.subplot(2,len(ind),ifig)
    plt.imshow(np.transpose(np.transpose(wc_out[ii,:,jy,:])), vmin=0.1, vmax=0.45, cmap='jet')
    plt.colorbar()
    plt.title('WC')
    ifig += 1

plt.show()
