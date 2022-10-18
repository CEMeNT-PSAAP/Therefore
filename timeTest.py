from tkinter import X
import numpy as np
#np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt
import therefore

#import mcdc
import numpy as np
import h5py

def t2p(time):
    return(int((time/max_time)*N_time))


# =============================================================================
# Therefore setup
# =============================================================================

data_type = np.float64

L = 10
dx = 0.1
N_mesh = int(L/dx)
xsec = 0.25
ratio = 0
scattering_xsec = xsec*ratio
source_mat = 0
N_angle = 4

dt = 2
max_time = 10

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
              'left_in_mag': 0.5,
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

'''
# =============================================================================
# MC/DC setup
# =============================================================================

# Set materials
m = mcdc.material(capture = np.array([xsec]),
                scatter = np.array([[0]]),
                fission = np.array([0.0]),
                nu_p    = np.array([0.0]),
                speed   = np.array([1.0]))

# Set surfaces
s1 = mcdc.surface('plane-x', x=-1E10, bc="vacuum")
s2 = mcdc.surface('plane-x', x=1E10,  bc="vacuum")

# Set cells
mcdc.cell([+s1, -s2], m)
mcdc.source(x=[0.0,1.0], time=np.array([0,1]), isotropic=True)

# Setting
mcdc.setting(N_particle=1E6)

# =============================================================================
# Running it
# =============================================================================

    
mcdc.tally(scores=['flux'], x=np.linspace(0, 1, 201), t=np.linspace(0, max_time, N_time+1))
mcdc.run()
'''

[sfMB, current, spec_rads] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'OCI_MB')
[sfEuler, current, spec_rads] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 1, 'OCI')

'''
with h5py.File('output.h5', 'r') as f:
    sfRef = f['tally/flux/mean'][:]
    t     = f['tally/grid/t'][:]
sfRef = np.transpose(sfRef)

for j in range(sfMB.shape[1]):
    sfEuler[:,j] = sfEuler[:,j]/max(sfEuler[:,j])
    sfMB[:,j] = sfMB[:,j]/max(sfMB[:,j])

for j in range(sfRef.shape[1]):
    sfRef[:,j] = sfRef[:,j]/max(sfRef[:,j])
'''
v=1
# analitic soultion
def psi(x,t):
    s = x + v*t
    return(0.5*np.exp((-xsec*s)/(1+angles_gq[1])))

x = np.linspace(0, L, int(N_mesh*2))
fig, axs = plt.subplots(N_time-1)
for i in range(1,N_time):
    axs[i-1].plot(x, sfEuler[:,i], label='euler')
    axs[i-1].plot(x, sfMB[:,i], label='mb')
    axs[i-1].plot(x, psi(x, i*dt), label='ref')
    axs[i-1].set_title('time {0}'.format(i*dt))
plt.tight_layout()

for ax in axs.flat:
    ax.label_outer()

plt.savefig('timeTest.png')

'''
x = np.linspace(0, L, int(N_mesh*2))
zippy = np.zeros(N_mesh*2)
plt.figure(4)
#plt.plot(x, ang_flux[:,0], label='MB 0')
#plt.plot(x, ang_flux[:,1], label='MB 1')
#plt.plot(x, ang_flux[:,2], label='MB 2')
#plt.plot(x, ang_flux[:,3], label='MB 3')
plt.plot(x, scalar_flux[:,0], label='Euler 0')
plt.plot(x, scalar_flux[:,1], label='Euler 1')
plt.plot(x, scalar_flux[:,2], label='Euler 2')
plt.plot(x, scalar_flux[:,3], label='Euler 3')
#plt.plot(x, zippy, '-r')
plt.xlabel('Distance [cm]')
plt.ylabel('Scalar Flux [units of scalar flux]')
plt.title('Scalar Flux through time')
plt.legend()
plt.show()'''
