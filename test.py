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

L = 4
dx = 1
xsec = 10
ratio = 0.9
scattering_xsec = xsec*ratio
source_mat = 0
source_a = 2
N_mesh = int(L/dx)

dx_mesh = dx*np.ones(N_mesh, data_type)
xsec_mesh = xsec*np.ones(N_mesh, data_type)
xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
source_mesh = source_mat*np.ones(N_mesh, data_type)

psi_in = source_mat / (xsec*(1-ratio)/2)
print(psi_in)


sim_perams = {'data_type': data_type,
              'N_angles': 2,
              'L': L,
              'N_mesh': N_mesh,
              'boundary_condition_left':  'incident_iso',
              'boundary_condition_right': 'incident_iso',
              'left_in_mag': 10,
              'right_in_mag': 10,
              'left_in_angle': .3,
              'right_in_angle': 1,
              'max loops': 10000}

[scalar_flux, current, spec_rad, conver] = therefore.OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
[scalar_flux2, current2, spec_rad2, conver2] = therefore.SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)

print('')
print('Did the implementaiton converge?'.format(conver))
print('Spectral radius of the lst run'.format(spec_rad))
print()
print()

print('>>>>OCI angular flux<<<<<<<<')
print(scalar_flux)
print()
print()
print('>>>>SI angular flux<<<<<<<<')
print(scalar_flux2)
print()
print()
print('>>>>OCI scalar flux<<<<<<<<')
print(np.sum(scalar_flux, axis=0))
print()
print()
print('>>>>SI scalar flux<<<<<<<<')
print(np.sum(scalar_flux2, axis=0))
print()
print()



f=1
X = np.linspace(0, L, int(N_mesh*2))
plt.figure(f)
<<<<<<< HEAD
plt.plot(X, scalar_flux[0,:], '-*k', label='OCI, 0')
plt.plot(X, scalar_flux[1,:], '--*k', label='OCI, 1')
plt.plot(X, scalar_flux2[0,:], '-r', label='SI, 0')
plt.plot(X, scalar_flux2[1,:], '--r', label='SI, 1')
plt.title('Test Flux')
plt.xlabel('Distance')
plt.ylabel('Angular Flux')
plt.legend()
plt.show()

f=2
#X = np.linspace(0, L, int(N_mesh*2))
X = np.arange(0, L, dx/2)
print(X)
plt.figure(f)
plt.plot(X, np.sum(scalar_flux, axis=0), '-*k', label='OCI')
plt.plot(X, np.sum(scalar_flux2, axis=0), '-r', label='SI')
plt.title('Test Flux')
plt.xlabel('Distance')
plt.ylabel('Scalar Flux')
plt.legend()
plt.show()


'''
=======
plt.plot(X, scalar_flux[0,:],  '-*k',  label='OCI 1')
plt.plot(X, scalar_flux[1,:],  '--*k', label='OCI 2')
plt.plot(X, scalar_flux2[0,:], '-r',  label='SI 1')
plt.plot(X, scalar_flux2[1,:], '--r', label='SI 2')
plt.title('Test Flux')
plt.xlabel('Distance')
plt.ylabel('Angular Flux')
plt.savefig('Test Angular flux')
'''
f+=1
X = np.linspace(0, L, int(N_mesh*2))
plt.figure(f)
plt.plot(X, scalar_flux[0,:], '-*k')
plt.plot(X, scalar_flux[1,:], '--*k')
plt.plot(X, scalar_flux2[0,:], '-r')
plt.plot(X, scalar_flux2[1,:], '--r')
plt.title('Test Flux')
plt.xlabel('Distance')
plt.ylabel('Scalar Flux')
plt.savefig('Test Scalar flux')
>>>>>>> 800494d79e530e7c85a7d478be7a974c1f69d86c
f+=1
plt.figure(f)
plt.title('Infinte Med')
plt.xlabel('Distance')
plt.ylabel('Current')
plt.ylim([-1,1])
flatLinePlot(X, current, '-k')
flatLinePlot(X, current2, '-r')
plt.show()
#launch source itterations #SourceItteration
'''
