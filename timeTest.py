import numpy as np
#np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt
import therefore

def flatLinePlot(x, y, dat):
    for i in range(y.size):
        xx = x[i:i+2]
        yy = [y[i], y[i]]
        plt.plot(xx, yy, dat)

data_type = np.float64

L = 1
dx = 0.01
N_mesh = int(L/dx)
xsec = 10
ratio = 0.0
scattering_xsec = xsec*ratio
source_mat = 0
N_angle = 2

dt = 0.1
max_time = .4

N_time = int(max_time/dt)
N_ans = 2*N_mesh

dx_mesh = dx*np.ones(N_mesh, data_type)
xsec_mesh = xsec*np.ones(N_mesh, data_type)
xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
source_mesh = source_mat*np.ones([N_mesh], data_type)

psi_in = source_mat / (xsec*(1-ratio)/2)
#print(psi_in)

[angles_gq, weights_gq] = np.polynomial.legendre.leggauss(N_angle)

#setup = np.linspace(0, np.pi, 2*N_mesh)
inital_scalar_flux = np.zeros(2*N_mesh)

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
              'max loops': 1000,
              'velocity': 1,
              'dt': dt,
              'max time': max_time,
              'N_time': N_time,
              'offset': 0,
              'ratio': ratio,
              'tolerance': 1e-9}

#inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source, backend='OCI_MB'
[scalar_flux, current, spec_rads] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 1, 'OCI')
[ang_flux, current, spec_rads] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'OCI_MB')

x = np.linspace(0, L, int(N_mesh*2))
plt.figure(4)
plt.plot(x, ang_flux[1,:,0], '-k', label='MB 0')
plt.plot(x, ang_flux[1,:,1], '--k', label='MB 1')
plt.plot(x, ang_flux[1,:,2], '-*k', label='MB 2')
plt.plot(x, ang_flux[1,:,3], '-^k', label='MB 3')
plt.plot(x, scalar_flux[1,:,0], '-b', label='Euler 0')
plt.plot(x, scalar_flux[1,:,1], '--b', label='Euler 1')
plt.plot(x, scalar_flux[1,:,2], '-*b', label='Euler 2')
plt.plot(x, scalar_flux[1,:,3], '-^b', label='Euler 3')
plt.xlabel('Distance [cm]')
plt.ylabel('Scalar Flux [units of scalar flux]')
plt.title('First time step of transient methods')
plt.legend()
plt.show()