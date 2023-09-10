"""
    Validate Serghei against Hydrus results reported in Beegum VZJ(2018)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

#   Output folder
fdir = ['out-omp4/','out-omp2-mpi2/']
lgd = ['OMP4','OMP2-MPI2','Hydrus']
#   Time index to plot
tt = [1]
#   Domain dimensions
L = 4000
H = 15
dz = 0.25
dim = [40, 1, 60]
N = int(dim[0]*dim[1]*dim[2])

"""
    --------------------------------------------------------------------
                          Read HYDRUS results
    --------------------------------------------------------------------
"""
hydrus = np.genfromtxt('HYDRUS2D.csv',delimiter=',',skip_header=1)

"""
    --------------------------------------------------------------------
                         Read Serghei results
    --------------------------------------------------------------------
"""
h_out = []
wc_out = []
for ff in range(len(fdir)):
    fname = fdir[ff] + 'output_subsurface.nc'
    fid = netcdf.NetCDFFile(fname,'r')
    h_out.append(fid.variables['hd'])
    wc_out.append(fid.variables['wc'])
    fid.close()
    print(' Shape of the output fields: ', np.shape(h_out[-1]))

wt_out = []
x = np.linspace(0, L, dim[0])
ioutput = 0
for ff in range(len(fdir)):
    for idx in range(len(tt)):
        serghei = np.transpose(np.transpose(h_out[ioutput][tt[idx],:,0,:]))
        wt = np.zeros((dim[0]))
        for ii in range(dim[0]):
            for kk in range(dim[2]-1):
                h1 = serghei[kk,ii]
                h2 = serghei[kk+1,ii]
                z1 = H - (kk+0.5)*dz
                z2 = H - (kk+1.5)*dz
                if h1 < 0.0 and h2 >= 0:
                    wt[ii] = z1 - h1*(z1-z2)/(h1-h2)
                    break
        wt_out.append(wt)
        ioutput += 1

"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""
fs = 10
co = ['b','g','r']
plt.figure(1, figsize=[8, 5])
ioutput = 0
for ff in range(len(fdir)):
    for idx in range(len(tt)):
        plt.plot(x, wt_out[ioutput], '-', color=co[ff])
        ioutput += 1
plt.scatter(hydrus[:,0], hydrus[:,1],  marker='o', facecolor='None', edgecolor='r')
plt.legend(lgd,fontsize=fs)
plt.xlim([0,L])
plt.xlabel('X [m]',fontsize=fs)
plt.ylabel('Water Table [m]',fontsize=fs)

# plt.savefig('watertable.png')

plt.show()
