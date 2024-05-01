"""
    Plot results for the tsunami test case

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file


fdir = ['out-sw/','out-pc/','out-mp/']
ind = [0,10,20,30]
jy = -1

h = []
z = []
qss = []
vol = []
vgw = []
vss = []
for ff in range(len(fdir)):
    fname = fdir[ff] + 'output.nc'
    fid = netcdf_file(fname,'r')
    h.append(fid.variables['h'])
    z.append(fid.variables['z'])
    if ff > 0:
        qss.append(fid.variables['qss'])
    else:
        qss.append([])
    fid.close()
    # get surface output
    fname = fdir[ff] + 'domainTimeSeries.out'
    vol.append(np.genfromtxt(fname, delimiter=' ', skip_header=1))
    # get subsurface output
    if ff > 0:
        fname = fdir[ff] + 'SubsurfaceTimeSeries.out'
        vgw.append(np.genfromtxt(fname, delimiter=' ', skip_header=1))
    else:
        vgw.append([])
    # get total sw-gw exchange volume [time, v_ss, v_tot]
    if ff > 0:
        tmp = []
        for ii in range(1,np.shape(vgw[-1])[0]):
            tmp.append([vgw[-1][ii,0], vgw[-1][ii,2]*(vgw[-1][ii,0]-vgw[-1][ii-1,0]), vgw[-1][ii,1]])
        vss.append(np.array(tmp))
    else:
        vss.append([])

print(' Shape of the output fields: ', np.shape(h))

"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""
cm2inch = 1.0/2.54
co = ['b','g','m','r']
ls = ['-','--',':']

plt.figure(1, figsize=[16*cm2inch, 18*cm2inch])
plt.subplot(2,1,1)
plt.plot(z[0][0,:],'k-')
for ff in range(len(fdir)):
    for ii in range(len(ind)):
        plt.plot(h[ff][ind[ii],0,:] + z[ff][0,:], color=co[ii], linestyle=ls[ff])
plt.title('Free surface elevation')
plt.xticks([],[])
plt.subplot(2,1,2)
for ff in range(len(fdir)):
    if ff > 0:
        plt.plot(abs(np.nansum(qss[ff][:,0,:],0)), color=co[ff], linestyle=ls[ff-1])
        # plt.scatter(np.linspace(0,100,100),np.nanmean(qss[ff][:,0,:],0), color=co[ii], linestyle=ls[ff-1])
    # plt.ylim([-1e-2, 1e-2])
    plt.yscale('log')
plt.xlabel('X [m]')
plt.ylabel('Flow Rate [m/s]')
plt.title('Mean SW-GW exchange rate [m/s]')


plt.figure(2, figsize=[16*cm2inch, 18*cm2inch])
plt.subplot(2,1,1)
plt.plot(vol[0][:,0], (vol[0][:,1]-vol[0][0,1])/vol[0][0,1], color=co[0])
plt.title('Surface water mass error')
plt.xlabel('Time [sec]')
plt.subplot(2,1,2)
for ff in range(len(fdir)):
    if ff > 0:
        plt.plot(vol[ff][:,0], ((vol[ff][:,1]+vgw[ff][:,1])-(vol[ff][0,1]+vgw[ff][0,1]))/(vol[ff][0,1]+vgw[ff][0,1]), color=co[ff])
plt.title('SW+GW water mass error')
plt.xlabel('Time [sec]')


plt.figure(3, figsize=[16*cm2inch, 18*cm2inch])

plt.subplot(2,1,1)
for ff in range(len(fdir)):
    if ff > 0:
        plt.plot(vss[ff][:,0], vss[ff][:,1], color=co[ff])
plt.subplot(2,1,2)
for ff in range(len(fdir)):
    if ff > 0:
        plt.plot(vol[ff][:,0], ((vol[ff][:,1]+vgw[ff][:,1])-(vol[ff][0,1]+vgw[ff][0,1])), color=co[ff])
    else:
        plt.plot(vol[ff][:,0], (vol[ff][:,1]-vol[ff][0,1]), color=co[ff])



# for ii in range(len(ind)):
#     plt.plot(qss[ind[ii],0,:], color=co[ii])

# for ii in range(len(ind)):
#     plt.subplot(2,2,ii+1)
#     plt.imshow(h[ind[ii],:,0,:], cmap='jet')


plt.show()
