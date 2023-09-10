"""
    Plot outflow for the sloping plane
"""

import numpy as np
import matplotlib.pyplot as plt

fdir = ['out-k5','out-k6','out-k7']
lgd = ['SERGHEI','ParFlow','CATHY']

"""
    Sulis2010 Results:
    (Abbrev. K5: Ks=1.17e-5 m/s, etc.)
    pflow = [[K5:t],[K5:Q],[K6:t],[K6:Q],[K7:t],[K7:Q]]

"""
pflow7 = [[35.5,63.6,84.0,106.4,143.3,199.2,209.6,222.4,244.9,285.9],
    [0.16,2.39,4.84,7.19,8.9,9.53,7.49,5.44,3.08,1.15]]
cathy7 = [[19.4,48.1,66.5,91.0,136.4,199.6,211.6,224.8,248.2,279.9],
    [0.09,2.51,5.02,7.62,9.34,9.73,7.57,5.47,3.05,1.47]]
pflow7 = np.array(pflow7)
cathy7 = np.array(cathy7)

pflow6 = [[61.5,102.6,132.9,160.0,198.8,212.7,237.3,275.4],
    [0.08,1.67,3.97,5.74,7.06,5.13,2.92,1.44]]
cathy6 = [[40.4,79.6,107.6,147.8,198.1,212.3,232.9,277.3],
    [0.08,1.74,3.82,6.04,7.21,5.15,2.71,0.69]]
pflow6 = np.array(pflow6)
cathy6 = np.array(cathy6)

pflow5 = [[122.4,149.1,167.4,197.3,211.5,228.9,253.6,291.5],
    [0.23,3.0,5.83,9.12,7.43,5.0,3.0,1.47]]
cathy5 = [[122.5,151.0,168.3,184.8,200.2,217.4,236.4,264.4],
    [0.05,2.39,4.96,7.35,8.77,6.34,4.25,2.38]]
pflow5 = np.array(pflow5)
cathy5 = np.array(cathy5)

data = []
for ff in range(len(fdir)):
    fname = fdir[ff]+'/domainTimeSeries.out'
    fid = open(fname, 'r')

    data.append(np.genfromtxt(fname, delimiter=' ', skip_header=True))

plt.figure(1, figsize=[14, 6])
ff = 0
plt.subplot(1,3,1)
plt.plot(data[ff][:,0]/60, data[ff][:,4]*60, color='b')
plt.scatter(pflow5[0,:], pflow5[1,:], marker='x', color='r')
plt.scatter(cathy5[0,:], cathy5[1,:], marker='^', facecolor='None', edgecolor='k')
plt.xlabel('Time [min]')
plt.ylabel('Flow Rate [m3/min]')
plt.ylim([0,10])
plt.title('Ks = 1.17e-5 m/s')
ff += 1

plt.subplot(1,3,2)
plt.plot(data[ff][:,0]/60, data[ff][:,4]*60, color='b')
plt.scatter(pflow6[0,:], pflow6[1,:], marker='x', color='r')
plt.scatter(cathy6[0,:], cathy6[1,:], marker='^', facecolor='None', edgecolor='k')
plt.xlabel('Time [min]')
plt.ylabel('Flow Rate [m3/min]')
plt.ylim([0,10])
plt.title('Ks = 1.17e-6 m/s')
plt.legend(lgd)
ff += 1

plt.subplot(1,3,3)
plt.plot(data[ff][:,0]/60, data[ff][:,4]*60, color='b')
plt.scatter(pflow7[0,:], pflow7[1,:], marker='x', color='r')
plt.scatter(cathy7[0,:], cathy7[1,:], marker='^', facecolor='None', edgecolor='k')
plt.xlabel('Time [min]')
plt.ylabel('Flow Rate [m3/min]')
plt.ylim([0,10])
plt.title('Ks = 1.17e-7 m/s')
ff += 1

plt.show()
