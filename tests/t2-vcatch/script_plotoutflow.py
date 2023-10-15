"""
    Plot results for the low triangle test problem
"""

import numpy as np
import matplotlib.pyplot as plt

fdir = ['out-rain','out-norain']
lgd = ['SERGHEI','CATHY','ParFlow','Cast3M']

# in h and m3/h
'''
	WITH RAIN
'''
pflow1 = [[1.95, 10.89, 19.33],
    [0, 1076.37, 1086.26]]
cathy1 = [[2.427, 3.156, 4.0, 5.568, 9.19, 19.8, 19.985, 20.675, 42.584, 72.582, 111.824],
    [7.545, 191.89, 513.61, 842.88, 1102.9, 1115.12, 745.13, 407.97, 65.15, 31.80, 20.30]]
pflow1 = np.array(pflow1)
cathy1 = np.array(cathy1)
'''
	NO RAIN
'''
pflow2 = [[15.86, 18.82, 24.74, 32.98, 48.34, 64.34, 79.19, 95.06, 116.23],
    [0.1, 1.53, 2.98, 3.77, 3.76, 3.31, 2.9, 2.49, 2.07]]
cathy2 = [[17.93, 21.4, 26.54, 38.41, 53.76, 69.76, 87.31, 106.16],
    [0.03, 1.35, 2.62, 3.52, 3.58, 3.43, 3.18, 2.95]]
c3m2 = [[17.68, 23.37, 32.93, 47.78, 63.64, 81.57, 99.63, 115.25],
    [0.04, 1.32, 2.26, 2.47, 2.28, 2.02, 1.77, 1.58]]
pflow2 = np.array(pflow2)
cathy2 = np.array(cathy2)
c3m2 = np.array(c3m2)

data = []
for ff in range(len(fdir)):
    fname = fdir[ff]+'/domainTimeSeries.out'
    fid = open(fname, 'r')
    data.append(np.genfromtxt(fname, delimiter=' ', skip_header=True))

plt.figure(1, figsize=[8, 4])
ff = 0

plt.subplot(1,2,1)
plt.plot(data[ff][:,0]/3600, 2.0*data[ff][:,4]*3600, color='r')
plt.scatter(cathy1[0,:], cathy1[1,:], marker='^', facecolor='None', edgecolor='k')
plt.scatter(pflow1[0,:], pflow1[1,:], marker='o', facecolor='None', edgecolor='k')
plt.xlabel('Time [h]')
plt.ylabel('Flow Rate [m3/h]')
plt.xlim([0,120])
plt.ylim([0,1200])
ff += 1

plt.subplot(1,2,2)
plt.plot(data[ff][:,0]/3600, 2.0*data[ff][:,4]*3600, color='r')
plt.scatter(cathy2[0,:], cathy2[1,:], marker='^', facecolor='None', edgecolor='k')
plt.scatter(pflow2[0,:], pflow2[1,:], marker='o', facecolor='None', edgecolor='k')
plt.scatter(c3m2[0,:], c3m2[1,:], marker='s', facecolor='None', edgecolor='k')
plt.xlabel('Time [h]')
plt.ylabel('Flow Rate [m3/h]')
plt.xlim([0,120])
plt.ylim([0,10])
plt.legend(lgd)

plt.savefig('outflow.png',format='png')

plt.show()
