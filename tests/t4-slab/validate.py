"""
    Validation of the Superslab case
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

fdir = 'out/'
lgd = ['HGS','CATHY','ParFlow','SERGHEI']
wcs = 0.1

'''
    Load reference data
'''
data = {}
keys = ['parflow','cathy','hgs']
fields = ['discharge','ponding','wcT3','wcT6']
for ii in range(len(fields)):
    data[fields[ii]] = {}
    for jj in range(len(keys)):
        fname = 'data/' + fields[ii] + '-' + keys[jj] + '.csv'
        data[fields[ii]][keys[jj]] = np.genfromtxt(fname, delimiter=',', skip_header=1)

'''
    Load model data
'''
out = np.genfromtxt(fdir + 'domainTimeSeries.out', skip_header=1)
t = out[:,0]
v = out[:,1]
q = out[:,4]

fid = netcdf.NetCDFFile(fdir + 'output_subsurface.nc','r')
wc = fid.variables['wc']
fid.close()
print(np.shape(wc))
wc3 = wc[6,:,0,60] / wcs
wc6 = wc[12,:,0,60] / wcs
z = np.linspace(-4.95,-0.05,50)


plt.figure(1, figsize=[14, 10])
mk = ['^','o','s']
co = ['g','r','b']
ff = 0
#   diacharge
plt.subplot(2,2,1)
for ii in range(len(keys)):
    plt.scatter(data['discharge'][keys[ii]][:,0], data['discharge'][keys[ii]][:,1], marker=mk[ii], facecolor='None', edgecolor=co[ii])
plt.plot(t/3600, q*3600, 'k-')
plt.xlabel('Time [h]')
plt.ylabel('Flow Rate [m3/h]')
ff += 1
#   surface ponding
plt.subplot(2,2,2)
for ii in range(len(keys)):
    plt.scatter(data['ponding'][keys[ii]][:,0], data['ponding'][keys[ii]][:,1], marker=mk[ii], facecolor='None', edgecolor=co[ii])
plt.plot(t/3600, v, 'k-')
plt.xlabel('Time [h]')
plt.ylabel('Ponding [m3]')
plt.legend(lgd)
ff += 1
#   water content at 3h
plt.subplot(2,2,3)
for ii in range(len(keys)):
    plt.scatter(data['wcT3'][keys[ii]][:,0], data['wcT3'][keys[ii]][:,1], marker=mk[ii], facecolor='None', edgecolor=co[ii])
plt.plot(np.flipud(wc3), z, 'k-')
plt.ylabel('Z [m]')
plt.xlabel('Saturation at 3h')
ff += 1
#   water content at 6h
plt.subplot(2,2,4)
for ii in range(len(keys)):
    plt.scatter(data['wcT6'][keys[ii]][:,0], data['wcT6'][keys[ii]][:,1], marker=mk[ii], facecolor='None', edgecolor=co[ii])
plt.plot(np.flipud(wc6), z, 'k-')
plt.ylabel('Z [m]')
plt.xlabel('Saturation at 6h')
ff += 1


plt.show()
