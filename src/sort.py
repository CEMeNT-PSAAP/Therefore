import numpy as np
import matplotlib.pyplot as plt
from mms import mms


N_angles = 2
N_cells = 100
N_groups = 2
N_time = 5



file_name_base = 'afluxUnsorted'
file_ext = '.csv'

[angles, weights] = np.polynomial.legendre.leggauss(N_angles)

# matrix order of the whole ass problem
SIZE_problem = N_cells*N_angles*N_groups*4
# size of the cell blocks in all groups and angle
SIZE_cellBlocks = N_angles*N_groups*4
# size of the group blocks in all angle within a cell
SIZE_groupBlocks = N_angles*4
# size of the angle blocks within a group and angle
SIZE_angleBlocks = 4

x = np.genfromtxt('x.csv', dtype=np.float64, delimiter=',', skip_header=1)
dx = x[1]

af_wp = np.zeros((N_time*2, N_groups, N_angles, 2*N_cells))

assert (int(af_wp.size/N_time) == SIZE_problem)

sf_wp = np.zeros((N_time*2, N_groups, 2*N_cells))

# schuky-duck the angular flux together
for t in range(N_time):
    # import csv file 
    file = file_name_base+str(t)+file_ext
    af_raw = np.genfromtxt(file, dtype=np.float64, delimiter=',', skip_header=2)
    af_raw = af_raw[:,0]

    if (af_raw.size != SIZE_problem):
        print(">>>ERROR<<<")
        print("Shape mismatch")
        print("af_raw shape: {0}".format(af_raw.size))
        print("SIZE_problem: {0}".format(SIZE_problem))
        #assert (af_raw.size == SIZE_problem)

    for i in range(N_cells):
        for g in range(N_groups):
            for n in range(N_angles):
                index_start = (i*(SIZE_cellBlocks) + g*(SIZE_groupBlocks) + 4*n)

                af_wp[t*2  ,g,n,2*i]   = af_raw[index_start]
                af_wp[t*2  ,g,n,2*i+1] = af_raw[index_start+1]
                af_wp[t*2+1,g,n,2*i]   = af_raw[index_start+2]
                af_wp[t*2+1,g,n,2*i+1] = af_raw[index_start+3]

                sf_wp[t*2  ,g,2*i]   += weights[n] * af_raw[index_start]
                sf_wp[t*2  ,g,2*i+1] += weights[n] * af_raw[index_start+1]
                sf_wp[t*2+1,g,2*i]   += weights[n] * af_raw[index_start+2]
                sf_wp[t*2+1,g,2*i+1] += weights[n] * af_raw[index_start+3]

#x = np.linspace(0, 1, N_cells*2)



# mms computation
sf_mms = np.zeros( (N_time*2, N_groups, N_cells) )
x_mms = np.linspace(0,1,N_cells)
dt = 1.0
manSol = mms(1,1,1,1,1,4)
for t in range(N_time):
    for i in range(N_cells):
        for g in range(N_groups):
            for n in range(N_angles):
                sf_mms[t*2,g,i]   += weights[n] * manSol.group2afCONT(angles[n],x_mms[i],dt*t)

print(sf_mms)
plt.figure()
plt.plot(x[:,0], sf_wp[7,0,:], label='computed')
plt.plot(x_mms, sf_mms[4,0,:], label='mms')
#plt.plot(x[:,0], sf_wp[2,1,:], label='g1 -- no source')
#plt.plot(x[:,0], sf_wp[5,0,:], label='g1 -- no source')
#plt.plot(x[:,0], sf_wp[5,1,:], label='g1 -- no source')
#plt.plot(x[:,0], sf_wp[7,0,:], label='g1 -- no source')
#plt.plot(x[:,0], sf_wp[7,1,:], label='g1 -- no source')
plt.xlabel('Distance')
plt.ylabel('Sc Fl')
plt.title('Single region -- trouble shoot time step=1')
plt.legend()
plt.show()