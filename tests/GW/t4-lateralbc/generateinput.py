"""
	Generate SERGHEI input files for the 3D rainfall-infiltration problem
"""
import numpy as np
import matplotlib.pyplot as plt

'''
	User should define grid resolutions and domain dimensions here
'''
#   Domain width
W = 100.0
#	Domain height
H = 15.0
#   Grid resolutions [m] in x, y, z directions
dx = 100.0
dy = 100.0
dz_base = 0.1
dz = []
r = 1.05

z = 0.0
kk = 1
dz.append(dz_base)
while z < H:
	dz.append(np.round(dz_base * r**kk, 3))
	kk += 1
	z += dz[-1]

dz = np.array(dz)
print(dz)

'''
	No user change below
'''
#   Domain length and height
L = 4000.0
H = 15.0
slope = 4.0 / L
dim = [int(L/dx), int(W/dy)]
nz = len(dz)
print('Dimension of the domain : ',dim[0], dim[1], nz)

#   Create DEM
dem = np.zeros((dim[1], dim[0]))
dem1d = []
for ii in range(dim[0]):
	dem1d.append(ii*dx*slope)
dem1d = np.array(dem1d)
for jj in range(dim[1]):
	dem[jj,:] = dem1d
dem = np.fliplr(dem)
dem1d = dem[0,:]
#   Save DEM output
hdstr = 'ncols '+str(dim[0])+'\nnrows '+str(dim[1])+'\nxllcorner 0.0\nyllcorner 0.0\ncellsize '+str(dx)+'\nnodata_value -9999'
np.savetxt('dem.input', dem, delimiter=' ', header=hdstr, fmt='%.3f', comments='')

#   Create Initial head
nynz = int(dim[1]*nz)
hic = np.zeros((nynz, dim[0]))
wtleft = np.amax(dem1d) - H + 7.0
wtright = np.amin(dem1d) - H + 0.9
wt = np.linspace(wtleft, wtright, dim[0])
irow = 0
for kk in range(nz):
	for jj in range(dim[1]):
		for ii in range(dim[0]):
			z = dem1d[ii] - np.nansum(dz[:kk]) - 0.5*dz[kk]
			#z = dem1d[ii] - (kk+0.5)*dz
			h = np.maximum(wt[ii] - z, -1.25)
			hic[irow,ii] = h
		irow += 1
#   Save Initial head file
icstr = 'ncols '+str(dim[0])+'\nnrows '+str(dim[1])+'\nnodata_value -9999'
np.savetxt('head.input', hic, delimiter=' ', header=icstr, fmt='%.3f', comments='')

#   Create Boundary head
hbcleft = np.reshape(hic[:,0], (nz, dim[1]))
hbcright = np.reshape(hic[:,-1], (nz, dim[1]))
bcstr = 'ndata 1 '+str(dim[1])+' '+str(nz)
np.savetxt('headleft.input', hbcleft, delimiter=' ', header=bcstr, fmt='%.3f', comments='')
np.savetxt('headright.input', hbcright, delimiter=' ', header=bcstr, fmt='%.3f', comments='')


