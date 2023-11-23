"""
    Plot results for the low triangle test problem
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


font_prop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf', size=10)  

fdir = ['output-rain','output-norain']
lgd = ['ParFlow','CATHY','HGS','Cast3M','SERGHEI']

# in h and m3/h
'''
	WITH RAIN
'''
pf1 = pd.read_csv('./ponding/PF-ponding-rain.csv', header=None)
cy1 = pd.read_csv('./ponding/CATHY-ponding-rain.csv', header=None)
hgs1 = pd.read_csv('./ponding/HGS-ponding-rain.csv', header=None)
cast1 = pd.read_csv('./ponding/cast3m-ponding-rain.csv', header=None)
'''
	NO RAIN
'''
pf2 = pd.read_csv('./ponding/PF-ponding-norain.csv', header=None)
cy2 = pd.read_csv('./ponding/CATHY-ponding-norain.csv', header=None)
hgs2 = pd.read_csv('./ponding/HGS-ponding-norain.csv', header=None)
cast2 = pd.read_csv('./ponding/cast3m-ponding-norain.csv', header=None)

p=[]
for ff in range(len(fdir)):
    fname = fdir[ff]+'/domainTimeSeries.out'
    out = np.genfromtxt(fname, skip_header =1)
    p.append(out[:,1])
    t=out[:,0]


plt.figure(1, figsize=[9, 4])

plt.subplot(1,2,1)

plt.scatter(pf1.iloc[:, 0], pf1.iloc[:, 1], s=14, marker='o', facecolor='None', edgecolor='k')
plt.scatter(cy1.iloc[:, 0], cy1.iloc[:, 1], s=14, marker='^', facecolor='None', edgecolor='k')
plt.scatter(hgs1.iloc[:, 0], hgs1.iloc[:, 1], s=14, marker='s', facecolor='None', edgecolor='k')
plt.scatter(cast1.iloc[:, 0], cast1.iloc[:, 1], s=14, marker='*', facecolor='None', edgecolor='k')
plt.plot(t/3600, 2*p[0], color='k')
plt.title("Scenario 1", fontproperties=font_prop)
plt.xlabel('Time [h]', fontproperties=font_prop)
plt.ylabel('ponded storage $m^{3}$', fontproperties=font_prop)
plt.xticks(fontproperties=font_prop)
plt.yticks(fontproperties=font_prop)
ff += 1

plt.subplot(1,2,2)
plt.scatter(pf2.iloc[:, 0], pf2.iloc[:, 1], s=14, marker='o', facecolor='None', edgecolor='k')
plt.scatter(cy2.iloc[:, 0], cy2.iloc[:, 1], s=14, marker='^', facecolor='None', edgecolor='k')
plt.scatter(hgs2.iloc[:, 0], hgs2.iloc[:, 1], s=14, marker='s', facecolor='None', edgecolor='k')
plt.scatter(cast2.iloc[:, 0], cast2.iloc[:, 1], s=14, marker='*', facecolor='None', edgecolor='k')
plt.plot(t/3600, 2*p[1], color='k')
plt.title("Scenario 2", fontproperties=font_prop)
plt.xlabel('Time [h]', fontproperties=font_prop)
plt.ylabel('ponded storage $m^{3}$', fontproperties=font_prop)
plt.legend(lgd,prop=font_prop)

plt.xticks(fontproperties=font_prop)
plt.yticks(fontproperties=font_prop)

plt.savefig('ponding2.png',format='png',bbox_inches='tight',dpi=600)

plt.show()
