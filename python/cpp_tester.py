import numpy as np
#np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt
import therefore
from timeit import default_timer as timer

def t2p(time):
    return(int((time/max_time)*N_time))


# =============================================================================
# Therefore setup
# =============================================================================

data_type = np.float64

L = 1
dx = 0.1
N_mesh = int(L/dx)
xsec = 1
ratio = 0
scattering_xsec = xsec*ratio
source_mat = 1
N_angle = 2

v = 1

dt = 0.5
max_time = .5

N_time = int(max_time/dt)
N_ans = 2*N_mesh
dx_mesh = dx*np.ones(N_mesh, data_type)
xsec_mesh = xsec*np.ones(N_mesh, data_type)
xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
source_mesh = source_mat*np.ones([N_mesh], data_type)

psi_in = source_mat / (xsec*(1-ratio)/2)
#print(psi_in)

[angles_gq, weights_gq] = np.polynomial.legendre.leggauss(N_angle)

IC = 0

inital_scalar_flux = IC*np.ones(2*N_mesh)

inital_angular_flux = np.zeros([N_angle, N_ans], data_type)
total_weight = sum(weights_gq)
for i in range(N_angle):
    inital_angular_flux[i, :] = inital_scalar_flux / total_weight

sim_perams = {'data_type': data_type,
              'N_angles': N_angle,
              'L': L,
              'N_mesh': N_mesh,
              'boundary_condition_left':  'vacuum',
              'boundary_condition_right': 'vacuum',
              'left_in_mag': 0,
              'right_in_mag': 0,
              'left_in_angle': 0,
              'right_in_angle': 0,
              'max loops': 10000,
              'velocity': v,
              'dt': dt,
              'max time': max_time,
              'N_time': N_time,
              'offset': 0,
              #'ratio': ratio,
              'tolerance': 1e-9,
              'print': True}


start = timer()
print('OCI MB SCB Single big gpu')
[sfMB, current, spec_rads, time] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'Big') #OCI_MB_GPU
end = timer()
print(end - start)
print(sfMB.shape)
x = np.linspace(0,1,sfMB.shape[0])

plt.figure()
plt.plot(x,sfMB[:,1])
plt.show()