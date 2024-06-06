
"""
    Validate Serghei against Hydrus results reported in Beegum VZJ(2018)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib.font_manager as fm

font_prop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf', size=18)  
fdir = ['output-t3/']
lgd = ['SERGHEI','Hydrus']
#simulation Length
T = 157680000
#outfrequence
F = 259200
#number of time steps
tt = int(T/F)

L = 4000
H = 15
dz = 0.25
dim = [40, 1, 60]
N = int(dim[0]*dim[1]*dim[2])

"""
    --------------------------------------------------------------------
                          Read HYDRUS results
    --------------------------------------------------------------------
"""
hydrus1 = np.genfromtxt('HYDRUS2D.csv',delimiter=',',skip_header=1)
hydrus2 = np.genfromtxt('hydrus-time.csv',delimiter=',',skip_header=1)

"""
    --------------------------------------------------------------------
                         Read Serghei results
    --------------------------------------------------------------------
"""
h_out = []
wc_out = []
for ff in range(len(fdir)):
    fname = fdir[ff] + 'output_subsurface.nc'
    fid = netcdf.NetCDFFile(fname,'r')
    h_out=fid.variables['hd']
    fid.close()
    print(' Shape of the output fields: ', np.shape(h_out))

wt_out = []
x = np.linspace(0, L, dim[0])
serghei = (h_out[-1,:,0,:])
wt = np.zeros((dim[0]))
for ii in range(dim[0]):
    for kk in range(dim[2]-1):
        h1 = serghei[kk,ii]
        h2 = serghei[kk+1,ii]
        z1 = H - (kk+0.5)*dz
        z2 = H - (kk+1.5)*dz
        if h1 < 0.0 and h2 >= 0:
            wt[ii] = z1 - h1*(z1-z2)/(h1-h2)
            break
wt_out.append(wt)

hcentre=[]
x_time=[]
for i in range(tt):
    aa=h_out[i,:,0,19]
    for kk in range(dim[2]-1):
        h1 = aa[kk]
        h2 = aa[kk+1]
        z1 = H - (kk+0.5)*dz
        z2 = H - (kk+1.5)*dz
        if h1 < 0.0 and h2 >= 0:
            bb = z1 - h1*(z1-z2)/(h1-h2)
            hcentre.append(bb) 
            break
for i in range(tt):
    a=i*259200/(3600*24)
    x_time.append(a)
"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""
color1=[0/255,150/255,136/255]
plt.figure(figsize=(16, 5)) 
for i in range(2):
    plt.subplot(1, 2, i+1)
    if i == 0:
        plt.plot(x, wt_out[0], '-', color='k')    
        plt.scatter(hydrus1[:,0], hydrus1[:,1],  marker='^', facecolor='None', edgecolor=color1)
        plt.legend(lgd)   
    else:
        xticks = range(0, 5 * 365+1, 365)
        xtick_labels = [str(i) for i in range(0, 5 * 365+1, 365)]
        plt.plot(x_time, hcentre, '-', color='k')  
        plt.scatter(hydrus2[:,0], hydrus2[:,1],  marker='^', facecolor='None', edgecolor=color1)  
    plt.xlabel('X [m]')
    plt.ylabel('Water Table [m]')



plt.savefig('both.png',format='png',bbox_inches='tight',dpi=600)
plt.show()

