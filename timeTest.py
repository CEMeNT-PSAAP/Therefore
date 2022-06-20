"""
Created on Sun May 15 20:34:03 2022
@author: jacksonmorgan
"""

import numpy as np
import matplotlib.pyplot as plt
import therefore


def flatLinePlot(x, y, dat):
    for i in range(y.size):
        xx = x[i:i+2]
        yy = [y[i], y[i]]
        plt.plot(xx, yy, dat)

data_type = np.float64

L = 10
dx = 1
xsec = 10
ratio = 0
scattering_xsec = xsec*ratio
source_mat = 1
#source_a = 2
N_mesh = int(L/dx)
N_angle = 2

dt = .01
max_time = .03
N_time = int(max_time/dt)
N_ans = 2*N_mesh

dx_mesh = dx*np.ones(N_mesh, data_type)
xsec_mesh = xsec*np.ones(N_mesh, data_type)
xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
source_mesh = source_mat*np.ones([N_mesh], data_type)

#psi_in = source_mat / (xsec*(1-ratio)/2)
#print(psi_in)

sim_perams = {'data_type': data_type,
              'N_angles': N_angle,
              'L': L,
              'N_mesh': N_mesh,
              'boundary_condition_left':  'vacuum',
              'boundary_condition_right': 'vacuum',
              'left_in_mag': 10,
              'right_in_mag': 10,
              'left_in_angle': .3,
              'right_in_angle': 0,
              'max loops': 10000,
              'velocity': 1,
              'dt': dt,
              'max time': max_time,
              'N_time': N_time}

theta = 1 #for discrete diamond
[scalar_flux, current, spec_rads] = therefore.TimeLoop(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, theta)


print(scalar_flux.shape)

#print('')
#print('Did the implementaiton converge?'.format(conver))
#print('Spectral radius of the lst run'.format(spec_rad))

#print(scalar_flux)

co = ['-k','-r','-b','-g','-p','-y']

f=1
X = np.linspace(0, L, int(N_mesh*2+1))
plt.figure(f)
for t in range(N_time):
    flatLinePlot(X, scalar_flux[:,t], co[t])
plt.title('Infinte Med')
plt.xlabel('Distance')
plt.ylabel('Scalar Flux')
plt.ylim([0,1.25*np.max(scalar_flux)])
plt.show()

'''
f+=1
plt.figure(f)
plt.title('Infinte Med')
plt.xlabel('Distance')
plt.ylabel('Current')
plt.ylim([-1,1])
flatLinePlot(X, current)
plt.show()
#launch source itterations #SourceItteration
'''
