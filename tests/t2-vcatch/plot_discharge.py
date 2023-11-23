"""
    Plot results for the discharge at the outlet
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

font_prop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf', size=14)  

fdir = ['output-norain','output-rain']
lgd = ['ParFlow','CATHY','HGS','Cast3M','SERGHEI']

'''
	WITH RAIN
'''
dataname='discharge/'
pf2 = pd.read_csv(dataname +'PF-rain.csv', header=None)
cy2 = pd.read_csv(dataname +'CATHY-rain.csv', header=None)
hgs2 = pd.read_csv(dataname +'HGS-rain.csv', header=None)

'''
	NO RAIN
'''
pf1 = pd.read_csv(dataname +'PF-norain.csv', header=None)
cy1 = pd.read_csv(dataname +'CATHY-norain.csv', header=None)
hgs1 = pd.read_csv(dataname +'HGS-norain.csv', header=None)
cast1 = pd.read_csv(dataname +'cast3m-norain.csv', header=None)

data = []
for ff in range(len(fdir)):
    fname = fdir[ff]+'/domainTimeSeries.out'
    fid = open(fname, 'r')
    data.append(np.genfromtxt(fname, delimiter=' ', skip_header=True))
print(' Shape of the output fields: ', len(data[0][:,4]))
plt.figure(1, figsize=[12, 4])
ff = 0

'''
	edgecolor
'''
color1 = (19/255, 103/255, 158/255)
color2 = (171/255,58/255,41/255)
color3 = (208/255, 127/255, 44/255)
color4 = (111/255, 109/255, 161/255)

plt.subplot(1,2,1)
plt.scatter(pf1.iloc[:, 0], pf1.iloc[:, 1], s=14, marker='o', facecolor='None', edgecolor=color1)


plt.scatter(cy1.iloc[:, 0], cy1.iloc[:, 1], s=14, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(hgs1.iloc[:, 0], hgs1.iloc[:, 1], s=14, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast1.iloc[:, 0], cast1.iloc[:, 1], s=14, marker='*', facecolor='None', edgecolor=color4)
plt.plot(data[ff][:,0]/3600, 2.0*data[ff][:,4]*3600,color='k')
plt.title("Scenario 1", fontproperties=font_prop)
plt.xlabel('Time [h]', fontproperties=font_prop)
plt.ylabel('Flow Rate [$m^{3}/h$]', fontproperties=font_prop)
plt.xlim([0,120])
plt.ylim([0,10])
plt.legend(lgd,prop=font_prop)
plt.xticks(fontproperties=font_prop)
plt.yticks(fontproperties=font_prop)
ff += 1

plt.subplot(1,2,2)
plt.scatter(pf2.iloc[:, 0], pf2.iloc[:, 1], s=14, marker='o', facecolor='None', edgecolor=color1)
plt.scatter(cy2.iloc[:, 0], cy2.iloc[:, 1], s=14, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(hgs2.iloc[:, 0], hgs2.iloc[:, 1], s=14, marker='s', facecolor='None', edgecolor=color3)
plt.plot(data[ff][:,0]/3600, 2.0*data[ff][:,4]*3600, color='k')
plt.title("Scenario 2", fontproperties=font_prop)
plt.xlabel('Time [h]', fontproperties=font_prop)
plt.ylabel('Flow Rate [$m^{3}/h$]', fontproperties=font_prop)
plt.xlim([0,120])
plt.ylim([0,1200])


plt.xticks(fontproperties=font_prop)
plt.yticks(fontproperties=font_prop)

plt.savefig('outflow.png',format='png',bbox_inches='tight',dpi=600)
plt.show()
