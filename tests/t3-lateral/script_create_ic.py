import numpy as np
import matplotlib.pyplot as plt


'''
    Initial condition
'''
nx = 40
nz = 60
dwt = [7.0, 0.9]
bottom = -15.0
thickH = 15.0

bath = np.linspace(3.9,0.0,nx)
dz = thickH / nz
wt = np.linspace(dwt[0],dwt[1],nx)
head = np.zeros((nz,nx))

for ii in range(nx):
    for kk in range(nz):
        head[kk,ii] = np.maximum(wt[ii] - (kk+0.5)*dz, -1.25)

head = np.flipud(head)

np.savetxt('head.input',head,delimiter=' ',fmt='%4.2f')

'''
    Boundary condition
'''
hbcxm = []
hbcxp = []
for kk in range(nz):
    hbcxm.append(dwt[0] - (kk+0.5)*dz)
    hbcxp.append(dwt[1] - (kk+0.5)*dz)
hbcxm = np.flipud(np.array(hbcxm))
hbcxp = np.flipud(np.array(hbcxp))
np.savetxt('hbcxm.input',hbcxm,delimiter=' ',fmt='%4.2f')
np.savetxt('hbcxp.input',hbcxp,delimiter=' ',fmt='%4.2f')

'''
    Rainfall
'''
days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
hours = np.array(days) * 24.0
q = [0.0004,0.0004,0.001,0.0017,0.0024,0.0019,0.0026,0.0021,0.0013,0.0011,0.0006,0.0004]
q = np.array(q) * 1e3 / 24.0
rain = []

year = 0
month = 0
t = 0.0
while True:
    rain.append([t, q[month]])
    t += hours[month]-0.1
    rain.append([t, q[month]])
    t += 0.1
    month += 1
    if month == 12:
        month = 0
        year += 1
    if year == 5:
        break

rain = np.round(rain,3)
# np.savetxt('rain',np.array(rain),delimiter=' ',fmt='%5.3f')


plt.figure(1)
plt.imshow(head)
plt.colorbar()
plt.show()
