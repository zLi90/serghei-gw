"""
    Validation of the Superslab case
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib.font_manager as fm
font_prop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf', size=18)  

fdir = ['output/','output-asy/']
lgd = ['ParFlow','CATHY','HGS','Cast3M','SERGHEI', 'SERGHEI-asy']
wcs = 0.1

'''
t3x0 means t=3h and x=0m
'''
dataname = 'data/'
pft3x0 = pd.read_csv(dataname +'PF-t3-x0.csv', header=None)
pft3x8 = pd.read_csv(dataname +'PF-t3-x8.csv', header=None)
pft3x40 = pd.read_csv(dataname +'PF-t3-x40.csv', header=None)
pft6x0 = pd.read_csv(dataname +'PF-t6-x0.csv', header=None)


chyt3x0 = pd.read_csv(dataname +'CATHY-t3-x0.csv', header=None)
chyt3x8 = pd.read_csv(dataname +'CATHY-t3-x8.csv', header=None)
chyt3x40= pd.read_csv(dataname +'CATHY-t3-x40.csv', header=None)
chyt6x0 = pd.read_csv(dataname +'CATHY-t6-x0.csv', header=None)
chyt6x8 = pd.read_csv(dataname +'CATHY-t6-x8.csv', header=None)
chyt6x40 = pd.read_csv(dataname +'CATHY-t6-x40.csv', header=None)

HGSt3x0 = pd.read_csv(dataname +'HGS-t3-x0.csv', header=None)
HGSt3x8 = pd.read_csv(dataname +'HGS-t3-x8.csv', header=None)
HGSt3x40 = pd.read_csv(dataname +'HGS-t3-x40.csv', header=None)
HGSt6x0 = pd.read_csv(dataname +'HGS-t6-x0.csv', header=None)
HGSt6x8 = pd.read_csv(dataname +'HGS-t6-x8.csv', header=None)
HGSt6x40 = pd.read_csv(dataname +'HGS-t6-x40.csv', header=None)

cast3mt3x0 = pd.read_csv(dataname +'cast3m-t3-x0.csv', header=None)
cast3mt3x8 = pd.read_csv(dataname +'cast3m-t3-x8.csv', header=None)
cast3mt3x40 = pd.read_csv(dataname +'cast3m-t3-x40.csv', header=None)
cast3mt6x0 = pd.read_csv(dataname +'cast3m-t6-x0.csv', header=None)
cast3mt6x8 = pd.read_csv(dataname +'cast3m-t6-x8.csv', header=None)
cast3mt6x40 = pd.read_csv(dataname +'cast3m-t6-x40.csv', header=None)

color1 = (19/255, 103/255, 158/255)
color2 = (171/255,58/255,41/255)
color3 = (208/255, 127/255, 44/255)
color4 = (111/255, 109/255, 161/255)

'''
    Load model data
'''

wc3x0= []
wc3x8= []
wc3x40=[]
wc6x0= []
wc6x8= []
wc6x40= []
t = []
v= []
q= []
for ff in range(len(fdir)):

    out = np.genfromtxt(fdir[ff] + 'domainTimeSeries.out', skip_header=1)
    t.append(out[:,0])
    v.append(out[:,1])
    q.append(out[:,4])

    fname = fdir[ff] + 'output_subsurface.nc'
    fid = netcdf.NetCDFFile(fname,'r')
    wc = fid.variables['wc']
    print(np.shape(wc)) 
    fid.close()
    wc3x0.append( wc[6,:,0,0] / wcs)
    wc3x8.append (wc[6,:,0,31] / wcs)
    wc3x40.append(wc[6,:,0,159] / wcs)
    wc6x0.append( wc[12,:,0,0] / wcs)
    wc6x8.append(wc[12,:,0,31] / wcs)
    wc6x40.append( wc[12,:,0,159] / wcs)



z = np.linspace(-4.95,-0.05,50)


plt.figure(1, figsize=[16, 12])
mk = ['^','o','s']
co = ['g','r','b']
ff = 0

plt.subplot(2,3,1)
plt.scatter(pft3x0.iloc[:,0],pft3x0.iloc[:,1]*(-1), s=30, marker='o', facecolor='None', edgecolor=color1)
plt.scatter(chyt3x0.iloc[:,0],chyt3x0.iloc[:,1]*(-1), s=30, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(HGSt3x0.iloc[:,0],HGSt3x0.iloc[:,1]*(-1),s=30, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mt3x0.iloc[:,0],cast3mt3x0.iloc[:,1]*(-1), s=30, marker='*', facecolor='None', edgecolor=color4)
plt.plot(np.flipud(wc3x0[0]), z, 'k-')
plt.plot(np.flipud(wc3x0[1]), z, 'r--')
plt.title("X=0m", fontproperties=font_prop)
plt.legend(lgd,prop=font_prop,loc="lower left")
plt.text(0.05, 0.95, "T=3h", transform=plt.gca().transAxes, fontsize=18, verticalalignment='top',fontproperties=font_prop)



plt.subplot(2,3,2)
plt.scatter(pft3x8.iloc[:,0],pft3x8.iloc[:,1]*(-1), s=30, marker='o', facecolor='None', edgecolor=color1)
plt.scatter(chyt3x8.iloc[:,0],chyt3x8.iloc[:,1]*(-1), s=30, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(HGSt3x8.iloc[:,0],HGSt3x8.iloc[:,1]*(-1),s=30, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mt3x8.iloc[:,0],cast3mt3x8.iloc[:,1]*(-1), s=30, marker='*', facecolor='None', edgecolor=color4)
plt.plot(np.flipud(wc3x8[0]), z, 'k-')
plt.plot(np.flipud(wc3x8[1]), z, 'r--')
plt.title("X=8m", fontproperties=font_prop)


plt.subplot(2,3,3)
plt.scatter(pft3x40.iloc[:,0],pft3x40.iloc[:,1]*(-1), s=30, marker='o', facecolor='None', edgecolor=color1)
plt.scatter(chyt3x40.iloc[:,0],chyt3x40.iloc[:,1]*(-1), s=30, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(HGSt3x40.iloc[:,0],HGSt3x40.iloc[:,1]*(-1),s=30, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mt3x40.iloc[:,0],cast3mt3x40.iloc[:,1]*(-1), s=30, marker='*', facecolor='None', edgecolor=color4)
plt.plot(np.flipud(wc3x40[0]), z, 'k-')
plt.plot(np.flipud(wc3x40[1]), z, 'r--')
plt.title("X=40m", fontproperties=font_prop)


plt.subplot(2,3,4)
plt.scatter(pft6x0.iloc[:,0],pft6x0.iloc[:,1]*(-1), s=30, marker='o', facecolor='None', edgecolor=color1)
plt.scatter(chyt6x0.iloc[:,0],chyt6x0.iloc[:,1]*(-1), s=30, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(HGSt6x0.iloc[:,0],HGSt6x0.iloc[:,1]*(-1),s=30, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mt6x0.iloc[:,0],cast3mt6x0.iloc[:,1]*(-1), s=30, marker='*', facecolor='None', edgecolor=color4)
plt.plot(np.flipud(wc6x0[0]), z, 'k-')
plt.plot(np.flipud(wc6x0[1]), z, 'r--')
plt.text(0.05, 0.95, "T=6h", transform=plt.gca().transAxes, fontsize=18, verticalalignment='top',fontproperties=font_prop)

plt.subplot(2,3,5)
plt.scatter(chyt6x8.iloc[:,0],chyt6x8.iloc[:,1]*(-1), s=30, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(HGSt6x8.iloc[:,0],HGSt6x8.iloc[:,1]*(-1),s=30, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mt6x8.iloc[:,0],cast3mt6x8.iloc[:,1]*(-1), s=30, marker='*', facecolor='None', edgecolor=color4)
plt.plot(np.flipud(wc6x8[0]), z, 'k-')
plt.plot(np.flipud(wc6x8[1]), z, 'r--')



plt.subplot(2,3,6)
plt.scatter(chyt6x40.iloc[:,0],chyt6x40.iloc[:,1]*(-1), s=30, marker='^', facecolor='None', edgecolor=color2)
plt.scatter(HGSt6x40.iloc[:,0],HGSt6x40.iloc[:,1]*(-1),s=30, marker='s', facecolor='None', edgecolor=color3)
plt.scatter(cast3mt6x40.iloc[:,0],cast3mt6x40.iloc[:,1]*(-1), s=30, marker='*', facecolor='None', edgecolor=color4)
plt.plot(np.flipud(wc6x40[0]), z, 'k-')
plt.plot(np.flipud(wc6x40[1]), z, 'r--')



for i in range(6):
    plt.subplot(2,3,i+1)
    plt.ylabel('Z [m]', fontproperties=font_prop)
    plt.xlabel('Saturation [-]', fontproperties=font_prop)
    plt.xticks([0.2, 0.4, 0.6, 0.8, 1.0],fontproperties=font_prop)
    plt.yticks(fontproperties=font_prop)
    plt.xlim(0, 1.0)
    plt.ylim(-5, 0)


plt.show()
plt.savefig('slab-saturation.png',format='png',bbox_inches='tight',dpi=600)
