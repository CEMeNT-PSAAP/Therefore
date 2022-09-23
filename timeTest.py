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

L = 1
dx = 0.25
N_mesh = int(L/dx)
xsec = 10
ratio = 0.5
scattering_xsec = xsec*ratio
source_mat = 0
N_angle = 2

dt = 0.25
max_time = 1

N_time = int(max_time/dt)
N_ans = 2*N_mesh

dx_mesh = dx*np.ones(N_mesh, data_type)
xsec_mesh = xsec*np.ones(N_mesh, data_type)
xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
source_mesh = source_mat*np.ones([N_mesh], data_type)

psi_in = source_mat / (xsec*(1-ratio)/2)
#print(psi_in)

#setup = np.linspace(0, np.pi, 2*N_mesh)
inital_angular_flux = np.zeros([N_angle, 2*N_mesh])
in_mid = np.ones(N_angle)

xm = np.linspace(-L/2,L/2, N_ans+1)
inital_scalar_flux = np.zeros(2*N_mesh)


[angles_gq, weights_gq] = np.polynomial.legendre.leggauss(N_angle)

assert(inital_scalar_flux.size == N_ans)

inital_angular_flux = np.zeros([N_angle, N_ans], data_type)
total_weight = sum(weights_gq)
for i in range(N_angle):
    inital_angular_flux[i, :] = inital_scalar_flux / total_weight

sim_perams = {'data_type': data_type,
              'N_angles': N_angle,
              'L': L,
              'N_mesh': N_mesh,
              'boundary_condition_left':  'incident_iso',
              'boundary_condition_right': 'vacuum',
              'left_in_mag': 10,
              'right_in_mag': 10,
              'left_in_angle': .3,
              'right_in_angle': 0,
              'max loops': 10000,
              'velocity': 1,
              'dt': dt,
              'max time': max_time,
              'N_time': N_time,
              'offset': 0,
              'ratio': ratio}


[scalar_flux, current, spec_rads] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 1, 'SI')
[scalar_flux2, current, spec_rads] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 1, 'OCI')

x = np.linspace(0, L, int(N_mesh*2))
plt.figure(4)
plt.plot(x, scalar_flux[:,0] , '-k', label='SI 0')
plt.plot(x, scalar_flux[:,1] , '-k', label='SI 1')
plt.plot(x, scalar_flux[:,2] , '-k', label='SI 2')
plt.plot(x, scalar_flux[:,3] , '-k', label='SI 3')
plt.plot(x, scalar_flux2[:,0], '-r', label='OCI 0')
plt.plot(x, scalar_flux2[:,1], '-r', label='OCI 1')
plt.plot(x, scalar_flux2[:,2], '-r', label='OCI 2')
plt.plot(x, scalar_flux2[:,3], '-r', label='OCI 3')
plt.xlabel('Distance [cm]')
plt.ylabel('Scalar Flux [units of scalar flux]')
plt.title('First time step of transient methods')
plt.legend()
plt.show()

'''
plt.figure(4)
plt.plot(x, scalar_flux[:,10] , '-k', label='SI')
plt.plot(x, scalar_flux2[:,10], '-r', label='OCI')
plt.plot(x, therefore.azurv1_spav(x_eval, ratio, inital_offset+dt*10),'--k', label='AVURV1')
plt.plot(x, inital_scalar_flux, '-b', label='Initial Condition')
plt.xlabel('Distance [cm]')
plt.ylabel('Scalar Flux [units of scalar flux]')
plt.title('First time step of transient methods')
plt.xlim(4.75,5.25)
plt.legend()
plt.savefig('SF_test.png')
'''