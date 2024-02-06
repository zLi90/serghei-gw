"""
    Verify serghei subsurface solver with Heterogeneous infiltration problem

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

fdir = ['out-pc/']
ind = [12]
jy = -1

ref400 = [
0.4415754218204367, 2.933358923111257,
0.4739169433928397, 2.464189817698022,
0.5569724031934337, 1.948641883813547,
0.7014208520650437, 1.7254228798340188,
1.038895686641265, 1.4953967340983398,
1.5459094967731146, 1.3935573142781312,
2.0465870127035672, 1.3550825007723024,
3.003016536258608, 1.3560404560025536,
3.8002530820632487, 1.4410219218144746,
4.2948989915472735, 1.7138848091260537,
4.477615393794263, 1.9706798384349602,
4.51966166655788, 2.416575841177304,
4.555375788363631, 2.803436847951056]
ref0 = [
1.0881349488838485, 2.1851481352297926,
1.5132033697700795, 2.1942940512019944,
2.0204794790884764, 2.2198039905904507,
2.6435280113858886, 2.224611935412572,
3.211901027218433, 2.2105659449863366,
3.7316773022665712, 2.1740560812738776,
3.9631271723344543, 2.0062676019048054,
3.512464023897719, 1.7817546026763578,
2.9267939224609187, 1.7477656822427021,
2.30330365666101, 1.7440292066475682,
1.7237681053368037, 1.7632350995996888,
1.282654177103391, 1.8151836707918303,
1.0150800041202963, 1.9801323708261096]

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

print(' Shape of the output fields: ', np.shape(h_out[0]))
# print(h_out[6,:,:])

"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""
ref0 = np.reshape(np.array(ref0), (13,2))
ref400 = np.reshape(np.array(ref400), (13,2))
color = ['k','b','r']
ls = ['-','--']
fs = 9
color = [[34,69,162],[143,38,126],[141,159,209],[209,140,184]]
for ii in range(len(color)):
    for jj in range(3):
        color[ii][jj] = color[ii][jj] / 255

cm2inch = 1.0/2.54

img = plt.figure(1, figsize=[14*cm2inch, 10*cm2inch])
ax = img.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0, pos1.y0+0.01, pos1.width*1.1, pos1.height*1.1]
ax.set_position(pos2)

xx = np.linspace(0,5,50)
zz = np.linspace(0,3,30)

wc2d = h_out[0][ind[0],:,0,:]
ct = ax.contour(xx,zz,np.flipud(wc2d), levels=[-400,0], cmap='RdBu')
ax.clabel(ct, fmt='%2.1f', colors='k', fontsize=14)

ax.scatter(ref0[:,0], ref0[:,1], s=30, marker='s', facecolors='none', edgecolors=color[0])
ax.scatter(ref400[:,0], ref400[:,1], s=30, marker='s', facecolors='none', edgecolors=color[1])

plt.xlabel('X [m]', fontsize=fs)
plt.ylabel('Z [m]', fontsize=fs)
# plt.legend(['SERGHEI','Reference'], fontsize=fs)

plt.savefig("t3-results.eps", format='eps')

plt.show()
