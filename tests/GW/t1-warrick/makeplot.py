"""
    Verify serghei subsurface solver with Warrick's analytical solution

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

analytical = np.array([[0.0825, 0.165, 0.2475], [74.58, 75.20, 77.12],
    [61.32, 62.13, 64.54], [39.02, 40.04, 42.80]])

fdir = ['out-pc-kk/','out-mp-kk/','out-pc-kkx5/','out-mp-kkx5/']
ind = [1, 2, 4]
jy = -1

h_out = []
wc_out = []
z_out = []
for ff in range(len(fdir)):
    fname = fdir[ff] + 'output_subsurface.nc'
    fid = netcdf_file(fname,'r')
    h_out.append(fid.variables['hd'])
    wc_out.append(fid.variables['wc'])
    z_out.append(fid.variables['z'])

fid.close()

print(' Shape of the output fields: ', np.shape(h_out))
# print(h_out[6,:,:])

"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""
color = [[141,159,209],[209,140,184],[34,69,162],[143,38,126]]
for ii in range(len(color)):
    for jj in range(3):
        color[ii][jj] = color[ii][jj] / 255
ls = ['-','-','--','--']
lw = 1.0
fs = 10
cm2inch = 1.0/2.54


img = plt.figure(1, figsize=[16*cm2inch, 14*cm2inch])

#   Plot the water content profiles
plt.subplot(1,2,1)
ax = img.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0-0.025, pos1.y0+0.01, pos1.width*1.1, pos1.height*1.1]
ax.set_position(pos2)

ifig = 1
for ii in range(len(ind)):
    for ff in range(len(fdir)):
        wc1d = wc_out[ff][ind[ii],:,0,0]
        plt.plot(wc1d, z_out[ff][:], color=color[ff], linestyle=ls[ff], linewidth=lw)
    plt.plot(analytical[0,:], -(1.0-analytical[ii+1,:]/100.0), linestyle='None',
        marker='x', markerfacecolor='None', color='k')
    ifig += 1
plt.xlabel('Water Content', fontsize=fs)
plt.ylabel('Z [m]', fontsize=fs)
plt.xticks([0,0.1,0.2,0.3],[0,0.1,0.2,0.3],fontsize=fs)
plt.yticks([0,-0.5,-1.0],[0,-0.5,-1.0],fontsize=fs)

plt.text(0.2,-0.28,'3.25h',fontsize=fs)
plt.text(0.2,-0.4,'6.5h',fontsize=fs)
plt.text(0.2,-0.63,'13h',fontsize=fs)
plt.text(0.01,0,'(a)',fontsize=fs)


#   Zoom in
plt.subplot(1,2,2)
ax = img.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0+0.025, pos1.y0+0.01, pos1.width*1.1, pos1.height*1.1]
ax.set_position(pos2)

ifig = 1
for ii in range(len(ind)):
    for ff in range(len(fdir)):
        wc1d = wc_out[ff][ind[ii],:,0,0]
        plt.plot(wc1d, z_out[ff][:], color=color[ff], linestyle=ls[ff], linewidth=lw)
    plt.plot(analytical[0,:], -(1.0-analytical[ii+1,:]/100.0), linestyle='None',
        marker='x', markerfacecolor='None', color='k')
    ifig += 1
plt.xlabel('Water Content', fontsize=fs)

plt.xlim([0.24,0.25])
plt.ylim([-0.36,-0.35])
plt.xticks([0.24,0.245,0.25],[0.24,0.245,0.25],fontsize=fs)
plt.yticks([-0.35,-0.355,-0.36],[-0.35,-0.355,-0.36],fontsize=fs)
plt.legend(['PC $\Delta t_{max}$=0.4s','MP $\Delta t_{max}$=0.4s','PC $\Delta t_{max}$=2s','MP $\Delta t_{max}$=2s','Analytical'],fontsize=fs)
plt.text(0.2403,-0.3504,'(b)',fontsize=fs)

plt.savefig("t1-results.eps", format='eps')
plt.show()
