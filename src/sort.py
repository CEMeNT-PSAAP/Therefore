import numpy as np
import matplotlib.pyplot as plt


N_angles = 16
N_cells = 20
N_groups = 2
N_time = 1

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

af_wp = np.zeros((N_time*2, N_groups, N_angles, 2*N_cells))

assert (af_wp.size == SIZE_problem)

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

x = np.linspace(0, 1, N_cells*2)

plt.figure()
plt.plot(x, sf_wp[1,1,:])
plt.show()