"""
    Validate Serghei against embankment problem
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

fdir = ['out-pc/']

ind = [0, 1, 3, 7, 10, 15]
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
    h_out.append(fid.variables['hd'])
    wc_out.append(fid.variables['wc'])
    fid.close()
    print(' Shape of the output fields: ', np.shape(h_out[ff]))

"""
    --------------------------------------------------------------------
                                Make plot
    --------------------------------------------------------------------
"""
reference = [
-4.531322864677335, 3.748773830159951,
-15.0628727191683, 3.6395462289290372,
-23.782391553475126, 3.615393054981359,
-34.09995820949513, 3.560336241643345,
-5.507899591434807, 3.393174263931425,
-13.621477519195711, 3.263969238159056,
-21.929667090941045, 3.2015999218991342,
-29.66017582321282, 3.135048840541316,
-1.2955221235177063, 3.1251201992488746,
-6.526342596153988, 2.6962225237579713,
-14.698522770718707, 2.5493126298816478,
-22.350010056406475, 2.4204087324467847,
-1.5286396737672092, 2.535681416302364,
-6.613244443814828, 2.1128790114923772,
-11.71171898331625, 1.9781285492580876,
-16.869619390135366, 1.8004511425001626,
-0.7340896695420369, 2.107607203063167,
-1.07773709275466, 1.4706877688727715,
-2.4982123316746154, 0.7350288340695775,
-1.3641659516267879, 0.19375855374063727
]
reference = np.array(np.reshape(reference, (20,2)))

color = [[223,223,234],[141,159,209],[209,140,184],[34,69,162],[143,38,126],[1,1,1]]
for ii in range(len(color)):
    for jj in range(3):
        color[ii][jj] = color[ii][jj] / 255
fs = 9

cm2inch = 1.0/2.54

img = plt.figure(1, figsize=[10*cm2inch, 14*cm2inch])
ax = img.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0, pos1.y0+0.01, pos1.width*1.1, pos1.height*1.1]
ax.set_position(pos2)

zz = np.linspace(4,0,40)
for ii in range(len(ind)):
    plt.plot(h_out[ff][ind[ii],:40,0,49]*9810/1000, zz, color=color[ii])
plt.scatter(reference[:,0],reference[:,1],s=40,c='k',marker='x')
plt.xlabel('Pressure [kPa]',fontsize=fs)
plt.ylabel('Z [m]',fontsize=fs)
plt.xticks([-40,-30,-20,-10,0],[-40,-30,-20,-10,0],fontsize=fs)
plt.yticks([0,1,2,3,4],[0,1,2,3,4],fontsize=fs)
plt.legend(['0 day','0.1 days','0.3 days','0.7 days','1 day','1.5 days','GFDM'], fontsize=fs)

# plt.savefig("t5-results.eps", format='eps')

plt.show()
