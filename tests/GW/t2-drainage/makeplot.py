"""
    Verify serghei subsurface solver with Free drainage test

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

fdir = ['out-pc/','out-mp/']
ind = [1, 4, 20, 99]
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

reference = [
0.16425893408572997, -0.39427695811338115,
0.1759420835121738, -1.1633316609812052,
0.1744452287855937, -1.9149458455362671,
0.18963217105749047, -2.7046960876317505,
0.20056584951282244, -3.4206490345254372,
0.19831494799621946, -0.4074901765366672,
0.20935061904675256, -1.159327879573385,
0.21477473056186774, -1.8895373953443295,
0.23000590574990445, -2.6901377952248815,
0.24076542211528065, -3.4424840952795526,
0.23967180139554448, -0.4320695475377816,
0.2553758112043876, -1.1381594111705917,
0.2611677944766922, -1.8929081385937003,
0.2784915639025152, -2.7125345584518112,
0.2790772687231271, -0.42880342996067533,
0.29968045531271303, -1.1394463597806936,
0.3099754457810949, -1.8973536728399703,
0.3184381820244707, -2.6953616810633014,
0.31985985684911783, -3.439031976507309
]
reference = np.reshape(np.array(reference),(19,2))

"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""
color = [[141,159,209],[143,38,126]]
# color = [[141,159,209],[209,140,184],[34,69,162],[143,38,126]]
for ii in range(len(color)):
    for jj in range(3):
        color[ii][jj] = color[ii][jj] / 255
ls = ['-','--','--','--']
lw = 1.0
fs = 9
cm2inch = 1.0/2.54

img = plt.figure(1, figsize=[10*cm2inch, 14*cm2inch])
ax = img.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0+0.035, pos1.y0+0.01, pos1.width, pos1.height*1.1]
ax.set_position(pos2)

ifig = 1
for ii in range(len(ind)):
    for ff in range(len(fdir)):
        wc1d = wc_out[ff][ind[ii],:,0,0]
        plt.plot(wc1d, z_out[ff][:], color=color[ff], linestyle=ls[ff])
        # print(h_out[ff][ind[ii],90:,0,0])
    plt.plot(reference[:,0], reference[:,1], linestyle='None', marker='x', markerfacecolor='None', color='k')
    ifig += 1

plt.text(0.155,-4,'100 days',fontsize=fs)
plt.text(0.205,-4,'20 days',fontsize=fs)
plt.text(0.26,-4,'4 days',fontsize=fs)
plt.text(0.31,-4,'1 day',fontsize=fs)

plt.xlabel('Water Content', fontsize=fs)
plt.ylabel('Z [m]', fontsize=fs)
# plt.xticks([0,0.1,0.2,0.3],[0,0.1,0.2,0.3],fontsize=fs)
# plt.yticks([0,-0.5,-1.0],[0,-0.5,-1.0],fontsize=fs)
plt.legend(['PC','MP','Experiment'], fontsize=fs)

# plt.savefig("t2-results.eps", format='eps')
plt.show()
