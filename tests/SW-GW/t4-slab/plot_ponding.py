import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import netcdf
import matplotlib.font_manager as fm

font_prop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf',size=10)  

fdir = ['output/','output-asy1/']
lgd = ['ParFlow','CATHY','HGS','Cast3M','SERGHEI', 'SERGHEI-asy']
wcs = 0.1
dataname = 'data/'

pfp = pd.read_csv(dataname+'PF-slab-ponding.csv', header=None)
cathyp = pd.read_csv(dataname+'CATHY-slab-ponding.csv', header=None)
hgsp = pd.read_csv(dataname+'HGS-slab-ponding.csv', header=None)
cast3mp = pd.read_csv(dataname+'cast3m-slab-ponding.csv', header=None)

pff = pd.read_csv(dataname+'PF-slab-outflow.csv', header=None)
cathyf = pd.read_csv(dataname+'CATHY-slab-outflow.csv', header=None)
hgsf = pd.read_csv(dataname+'HGS-slab-outflow.csv', header=None)
cast3mf = pd.read_csv(dataname+'cast3m-slab-outflow.csv', header=None)

color1 = (19/255, 103/255, 158/255)
color2 = (171/255,58/255,41/255)
color3 = (208/255, 127/255, 44/255)
color4 = (111/255, 109/255, 161/255)

t = []
p= []
q= []
for ff in range(len(fdir)):

    out = np.genfromtxt(fdir[ff] + 'domainTimeSeries.out', skip_header=1)
    t.append(out[:,0])
    p.append(out[:,1])
    q.append(out[:,4])


plt.figure(1, figsize=[9, 4])
plt.rc('font', size=10)
ss=14
plt.subplot(1,2,1)
plt.scatter(pfp.iloc[:,0],pfp.iloc[:,1], s=ss, marker='o', facecolor='None', edgecolor=color1)
plt.scatter(cathyp.iloc[:,0],cathyp.iloc[:,1], s=ss, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(hgsp.iloc[:,0],hgsp.iloc[:,1],s=ss, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mp.iloc[:,0],cast3mp.iloc[:,1], s=ss, marker='*', facecolor='None', edgecolor=color4)
plt.plot(t[0]/3600, p[0], 'k-')
plt.plot(t[0]/3600, p[1], 'r--')
plt.xlabel('Time [h]',fontproperties=font_prop)
plt.ylabel('Ponding storage [$m^{3}$]',fontproperties=font_prop)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(fontproperties=font_prop)
plt.yticks(fontproperties=font_prop)
plt.legend(lgd,prop=font_prop,loc="upper right")


plt.subplot(1,2,2)
plt.scatter(pff.iloc[:,0],pff.iloc[:,1], s=ss, marker='o', facecolor='None', edgecolor=color1)
plt.scatter(cathyf.iloc[:,0],cathyf.iloc[:,1], s=ss, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(hgsf.iloc[:,0],hgsf.iloc[:,1],s=ss, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mf.iloc[:,0],cast3mf.iloc[:,1], s=ss, marker='*', facecolor='None', edgecolor=color4)
plt.plot(t[0]/3600, q[0]*3600, 'k-')
plt.plot(t[0]/3600, q[1]*3600, 'r--')
plt.xlabel('Time [h]',fontproperties=font_prop)
plt.ylabel('Flow Rate[$m^{3}$/h]',fontproperties=font_prop)
plt.xticks(fontproperties=font_prop)
plt.yticks(fontproperties=font_prop)



plt.savefig('ponding-flow.png',format='png',bbox_inches='tight',dpi=600)

plt.show()