import numpy as np
import matplotlib.pyplot as plt


N_angles = 2
N_cells = 170
N_groups = 1
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


dx_mesh = np.empty(N_cells)
N_region = np.array((8, 8, 4, 50, 100), int)
dx = np.array((0.25, 0.25, 0.25, 0.02, 0.02))
for i in range(N_region.size):
    LB = sum(N_region[:i])
    RB = sum(N_region[:i+1])
    dx_mesh[LB:RB] = dx[i]

x = np.zeros(N_cells*2)
for i in range(N_cells):
    x[2*i] = sum(dx_mesh[:i])
    x[2*i+1] = sum(dx_mesh[:i+1])

af_wp = np.zeros((N_time*2, N_groups, N_angles, 2*N_cells))

assert (int(af_wp.size/N_time) == SIZE_problem)

sf_wp = np.zeros((N_time*2, N_groups, 2*N_cells))

# schuky duck the angular flux together
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
        assert (af_raw.size == SIZE_problem)

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

plt.figure()
plt.plot(x, sf_wp[1,0,:], label='1')
plt.plot(x, sf_wp[3,0,:], label='2')
plt.plot(x, sf_wp[5,0,:], label='3')
plt.plot(x, sf_wp[6,0,:], label='4')
plt.xlabel('Distance')
plt.ylabel('Sc Fl')
plt.title('Trans Reeds -- trouble shoot')
plt.legend()
plt.show()