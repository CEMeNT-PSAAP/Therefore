from tkinter import X
import numpy as np
#np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt
import therefore
from timeit import default_timer as timer

#import mcdc
import numpy as np
#import h5py

def t2p(time):
    return(int((time/max_time)*N_time))


# =============================================================================
# Therefore setup
# =============================================================================

data_type = np.float64

L = 10
dx = 0.01
N_mesh = int(L/dx)
xsec = 0.25
ratio = 0.75
scattering_xsec = xsec*ratio
source_mat = 0
N_angle = 32

v = 1

BCl = 0.5

dt = 0.1
max_time = 5

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
              'left_in_mag': BCl,
              'right_in_mag': 10,
              'left_in_angle': .3,
              'right_in_angle': 0,
              'max loops': 10000,
              'velocity': v,
              'dt': dt,
              'max time': max_time,
              'N_time': N_time,
              'offset': 0,
              'ratio': ratio,
              'tolerance': 1e-9,
              'print': False}

start = timer()
print('OCI MB SCB Single big gpu')
[sfMB, current, spec_rads] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'BigGirl') #OCI_MB_GPU
end = timer()
print(end - start)

'''
start = timer()
print('OCI MB SCB Small GPU')
[sfMB, current, spec_rads] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'OCI_MB_GPU')
end = timer()
print(end - start)

'''

start = timer()
print('OCI MB SCB CPU')
[sfMB_trad, current, spec_rads] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'OCI_MB')
end = timer()
print(end - start)

'''
start = timer()
print('SI MB SCB')
[sfMBSi, current, spec_rads] = therefore.multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'SI_MB')
end = timer()
print(end - start)

'''

start = timer()
print('SI BE SCB')
[sfEuler, current, spec_rads, loops] = therefore.euler(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, 'SI')
end = timer()
print(end - start)


x = np.linspace(0, L, int(N_mesh*2))

fig,ax = plt.subplots()
    
ax.grid()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\phi$')
ax.set_title('Scalar Flux (ϕ)')

import matplotlib.animation as animation

line1, = ax.plot(x, sfMB[:,0], '-k',label="MB-SCB-OCI-GPU")
line2, = ax.plot(x, sfEuler[:,0], '-r',label="BE-SCB-SI")
line3, = ax.plot(x, sfMB_trad[:,0], '-g',label="MB-SCB-OCI-CPU")
text   = ax.text(8.0,0.75,'') 
ax.legend()
plt.ylim(-0.2, 1.5)

def animate(k):
    line1.set_ydata(sfMB[:,k])
    line2.set_ydata(sfEuler[:,k])
    line3.set_ydata(sfMB_trad[:,k])

    #line3.set_ydata(sfMBSi[:,k])
    text.set_text(r'$t \in [%.1f,%.1f]$ s'%(dt*k,dt*(k+1)))
    return line1#, line2

simulation = animation.FuncAnimation(fig, animate, frames=N_time)

writervideo = animation.PillowWriter(fps=250)
simulation.save('test_mb.gif') #saveit!