import numpy as np
#np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt
import therefore

import mcdc
import numpy as np
import h5py


# =============================================================================
# Therefore setup
# =============================================================================

data_type = np.float64

L = 1
dx = 0.01
N_mesh = int(L/dx)
xsec = 1
ratio = 0.0
scattering_xsec = xsec*ratio
source_mat = 0
N_angle = 2

dt = np.array([.1])
max_time = 1

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
              'dt': dt[0],
              'max time': max_time,
              'N_time': N_time,
              'offset': 0,
              'ratio': ratio,
              'tolerance': 1e-9}


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

mcdc.source(point=[1E-9,0.0,0.0], time=np.array([0,1]), isotropic=True)

# Tally
mcdc.tally(scores=['flux', 'flux-t'], x=np.linspace(0, 1, 201), t=np.linspace(0, max_time, N_time+1)) #np.arange(0, 0.4, dt)

# Setting
mcdc.setting(N_particle=1E5)

# =============================================================================
# Running it
# =============================================================================


errorMB = np.zeros(dt.size)
errorEuler = np.zeros(dt.size)

for i in range(dt.size):
    sim_perams['dt'] = dt[i]
    N_time = int(max_time/dt[i])
    #mcdc.tally(scores=['flux'], x=np.linspace(0, 1, 201), t=np.linspace(0, max_time, N_time+1))
    mcdc.run()

    [sfEuler, current, spec_rads] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 1, 'OCI')
    [sfMB, current, spec_rads] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'OCI_MB')

    with h5py.File('output.h5', 'r') as f:
        sfRef = f['tally/flux/mean'][:]
        t     = f['tally/grid/t'][:]
    
    print(t)
    sfRef = np.transpose(sfRef)

    for j in range(N_time):
        sfEuler[:,j] = sfEuler[:,j]/max(sfEuler[:,j])
        sfMB[:,j] = sfMB[:,j]/max(sfMB[:,j])
    for j in range(3):
        sfRef[:,j] = sfRef[:,j]/max(sfRef[:,j])
        
        

    print()
    print('dt slice')
    print(' -Solution Euler Shape:     {0}'.format(sfEuler.shape))
    print(' -Solution TMDB Shape:      {0}'.format(sfMB.shape))
    print(' -Reference Solution Shape: {0}'.format(sfRef.shape))

    assert (sfEuler.shape == sfMB.shape)
    #assert (sfRef.shape == sfEuler.shape)
    '''
    for j in range(N_time):
        errorMB[i] += np.linalg.norm(sfMB[:,j]-sfRef[:,j], ord == 2)
        errorEuler[i] += np.linalg.norm(sfEuler[:,j]-sfRef[:,j], ord == 2)
    errorMB[i] /= N_time
    errorEuler[i] /= N_time

    print(' -multi balance error:      {0}'.format(errorMB))
    print(' -euler error:              {0}'.format(errorEuler))
    '''


print()
x = np.linspace(0, L, int(N_mesh*2))

fig, axs = plt.subplots(4)
axs[0].plot(x, sfEuler[:,1], label='euler')
axs[0].plot(x, sfMB[:,1], label='mb')
axs[0].plot(x, sfRef[:,1], label='ref')
axs[0].plot(x, sfRef[:,0], label='ref')
axs[0].set_title('0.1 [s]')

axs[1].plot(x, sfEuler[:,2], label='euler')
axs[1].plot(x, sfMB[:,2], label='mb')
axs[1].plot(x, sfRef[:,2], label='ref')
axs[1].set_title('0.2 [s]')

axs[2].plot(x, sfEuler[:,3], label='euler')
axs[2].plot(x, sfMB[:,3], label='mb')
axs[2].plot(x, sfRef[:,3], label='ref')
axs[2].set_title('0.3 [s]')

print(' -Reference Solution Shape 2: {0}'.format(sfRef.shape))
axs[3].plot(x, sfRef[:,0], label='1')
axs[3].plot(x, sfRef[:,1], label='2')
axs[3].plot(x, sfRef[:,2], label='3')
axs[3].plot(x, sfRef[:,3], label='3')
axs[3].plot(x, sfRef[:,4], label='3')
axs[3].set_title('0.4 [s]')

for ax in axs.flat:
    ax.label_outer()

plt.show()

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
