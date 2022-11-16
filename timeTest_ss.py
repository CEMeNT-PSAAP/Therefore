from tkinter import X
import numpy as np
#np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt
import therefore

#import mcdc
import numpy as np
#import h5py

import scipy.special as sc

def t2p(time):
    return(int((time/max_time)*N_time))


def phi_(x,t):
    if x > v*t:
        return 0.0
    else:
        return 1.0/bound_mag * (xsec*x*(sc.exp1(xsec*v*t) - sc.exp1(xsec*x)) + \
                        np.e**(-xsec*x) - x/(v*t)*np.e**(-xsec*v*t))

def analitical(x, t):
    y = np.zeros(x.shape)
    for i in range(x.size):
        y[i] = phi_(x[i],t)
    return y

# =============================================================================
# Therefore setup
# =============================================================================

data_type = np.float64

L = 10
dx = 0.001
N_mesh = int(L/dx)
xsec = 0.25
ratio = 0.75
scattering_xsec = xsec*ratio
source_mat = 0
N_angle = 2
bound_mag = 1
v=1

dt = 2
max_time = 12

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
              'boundary_condition_right': 'reflecting',
              'left_in_mag': bound_mag,
              'right_in_mag': 10,
              'left_in_angle': .3,
              'right_in_angle': 0,
              'max loops': 10000,
              'velocity': 1,
              'dt': dt,
              'max time': max_time,
              'N_time': N_time,
              'offset': 0,
              'ratio': ratio,
              'tolerance': 1e-9,
              'print': True}

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
#
dt_t = np.array([ .1]) #[2, 1.5, 1, .75, .5, .25,
N_dt = dt_t.size

errorMB = np.zeros_like(dt_t)
errorEuler = np.zeros_like(dt_t)

def printer(i):
    print('')
    print('========================================')
    print('         cycle {0}'.format(i))
    print('========================================')
    print('')


#fig, axs = plt.subplots(dt_t.size)

for i in range(dt_t.size):

    printer(i)

    N_time = int(max_time/dt_t[i])

    sim_perams['dt'] = float(dt_t[i])
    sim_perams['N_time'] = N_time

    x = np.linspace(0, L, int(N_mesh*2))
    sfRef = analitical(x, max_time)
    sfRef[0] = 1

    [sfMB, current, spec_rads] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'SI')
    [sfEuler, current, spec_rads] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 1, 'OCI')

    errorMB[i] =    np.linalg.norm((sfMB[:,-1]-sfRef) / sfRef)
    errorEuler[i] = np.linalg.norm((sfEuler[:,-1]-sfRef) / sfRef)

    #axs[i].plot(x, sfEuler[:,-1], label='euler')
    #axs[i].plot(x, sfMB[:,-1], label='mb')
    #axs[i].plot(x, analitical(x, max_time), label='ref')
    #axs[i].set_title('dt {0}[s]'.format(dt_t[i]))
#plt.tight_layout()
#axs[-1].legend()
#plt.savefig('timeTest_overdt.pdf', dpi=200)

print(errorMB)
print(errorEuler)


plt.figure()
plt.plot(x, sfEuler[:,-1], label='euler')
plt.plot(x, sfMB[:,-1], label='mb')
#plt.plot(x, analitical(x, max_time), label='ref')
plt.title('time {0}[s]'.format(max_time))

plt.savefig('timeTest.pdf', dpi=200)


plt.figure()
plt.plot(dt_t, errorMB, label='MB')
plt.plot(dt_t, errorEuler, label='E')
plt.legend()
plt.savefig('order_est.pdf', dpi=200)





#[ss_sfMB, current, spec_rads, tf] = therefore.OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)


'''
ss_xRef = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10])
ss_sfRef = np.array([0.299999998, 0.268779374, 0.241460362, 0.217481669, 0.196369686, 0.177724223, 0.161206604, 0.146529743, 0.133449844, 0.121759478, 0.111281796, 0.101865694, 0.093381763, 0.085718903, 0.078781493, 0.072487008, 0.066764032, 0.061550582, 0.056792705, 0.052443293, 0.048461094, 0.044809867, 0.041457677, 0.038376299, 0.035540704, 0.032928639, 0.030520253, 0.028297791, 0.02624533, 0.024348547, 0.022594531, 0.020971611, 0.019469216, 0.018077747, 0.016788474, 0.015593439, 0.014485377, 0.013457643, 0.01250415, 0.011619316, 0.010798014])
'''



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


v=1

fig, axs = plt.subplots(int((N_time)/2))
for i in range(1, N_time, 2):
    it = int((i-1)/2)
    axs[it].plot(x, sfEuler[:,i], label='euler')
    axs[it].plot(x, sfMB[:,i], label='mb')
    axs[it].plot(x, analitical(x, i*dt), label='ref')
    axs[it].set_title('time {0}[s]'.format(i*dt))

plt.tight_layout()
axs[-1].legend()

for ax in axs.flat:
    ax.label_outer()

plt.savefig('timeTest.pdf', dpi=200)


x = np.linspace(0, L, int(N_mesh*2))
plt.figure(1)
#plt.plot(x, sfEuler[:,-1], label='euler, t=20')
#plt.plot(x, sfMB[:,-1], label='mb, t=20')
plt.plot(ss_xRef, ss_sfRef, label='SS excel')
plt.plot(x, ss_sfMB, label='SS OCI')
plt.title('Steady State (t=20[s]) v=1, sig=0.25, N=100, c=0')
plt.legend()

plt.savefig('timeTest_ss.png')


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
