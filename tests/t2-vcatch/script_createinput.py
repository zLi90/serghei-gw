"""
    Create inputs for the V-Catchment test problem
"""
import numpy as np
import matplotlib.pyplot as plt

'''
    Domain
'''
L = 100.0
W = 50.0
H = 5.0
delta = [1.0, 1.0, 0.2]
nx = int(L/delta[0])
ny = int(W/delta[1])
nz = int(H/delta[2])
wt = 2.0

savedem = False
saveRoughness = True
saveWaterTable = False

bath = np.zeros((ny, nx))
for ii in range(ny):
    bath[ii,:] = 0.05 * ii * delta[1]
    for jj in range(nx):
        bath[ii,jj] += 0.02 * jj * delta[0]

print(np.nanmax(bath), np.nanmin(bath))

bath2 = np.zeros((ny+5, nx))
bath2[5:,:] = bath + 0.05
for ii in range(5):
    bath2[ii,:] = bath[0,:]
if savedem:
    np.savetxt('dem.input',bath2,delimiter=' ',fmt='%.3f')

roughness = 0.626 * np.ones((ny,nx))
roughness2 = 6.26 * np.ones((ny+5, nx))
roughness2[5:,:] = roughness
if saveRoughness:
    np.savetxt('roughness.input',roughness2,delimiter=' ',fmt='%.3f')

wt = bath2 - 2.0
if saveWaterTable:
    np.savetxt('wt.input',wt,delimiter=' ',fmt='%.3f')


plt.imshow(bath2)
plt.show()
