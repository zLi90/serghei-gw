"""
	Compare simulation against measured data
"""
import numpy as np
import matplotlib.pyplot as plt

nsite = 6
ndays = 30
'''
	Load simulation results
'''
fout = 'out-sw/gauges.out'
gauges = np.genfromtxt(fout, delimiter='\t', skip_header=1)

'''
	Load validation data
'''
dt = 15 * 60
ndata = int(ndays*86400.0/dt)
fvad = 'data/validation/'
valid = []
for ii in range(1,nsite+1):
	tmp = np.genfromtxt(fvad+'waterlevel_0515-0731_nueces'+str(ii), delimiter=' ')
	tmp = tmp[:ndata]
	for jj in range(len(tmp)):
		if tmp[jj] < -100:
			tmp[jj] = np.nan
	valid.append(tmp)
tt = np.arange(0, dt*ndata, dt)

'''
	Make plots
'''
plt.figure(1)
ifig = 1
for ii in range(3):
	for jj in range(2):
		plt.subplot(3,2,ifig)
		plt.plot(tt/86400, valid[ifig-1], 'k-')
		plt.plot(gauges[:,0]/86400.0, gauges[:,int((ifig-1)*4+1)] + gauges[:,int(ifig*4)])
		plt.xlim([0,ndays])
		ifig += 1

plt.show()
