"""
Created on Sun May 15 20:34:03 2022
@author: jacksonmorgan
"""

import numpy as np
import matplotlib.pyplot as plt
import therefore

def analiticalSoultion(t, xsec, inital_flux, source, vel):
    angular_flux = inital_flux*np.exp(-xsec*vel*t) + (source/xsec)*(1-np.exp(-xsec*vel*t))
    
    return(angular_flux)
    
#\Psi (t) = \Psi(0) e^(-\sigma*v*t) + Q/(\sigma)*(1-np.exp**(-\sigma*v*t))

def flatLinePlot(x, y, dat):
    for i in range(y.size):
        xx = x[i:i+2]
        yy = [y[i], y[i]]
        plt.plot(xx, yy, dat)

data_type = np.float64

L = 10
dx = .01
xsec = 1
ratio = 0.9
scattering_xsec = xsec*ratio
source_mat = 0
#source_a = 2
N_mesh = int(L/dx)
N_angle = 2

dt = .1
max_time = 10
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
for i in range(int(.45*N_mesh*2), int(.55*N_mesh*2), 1):
    inital_angular_flux[:,i] = in_mid

#inital_angular_flux = np.array([[np.sin(setup)],[np.sin(setup)]]).reshape(2,200) #[:,:N_mesh]

#in_an_flux = 1
#inital_angular_flux = in_an_flux * np.ones([N_angle, 2*N_mesh])

sim_perams = {'data_type': data_type,
              'N_angles': N_angle,
              'L': L,
              'N_mesh': N_mesh,
              'boundary_condition_left':  'reflecting',
              'boundary_condition_right': 'reflecting',
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
[scalar_flux, current, spec_rads] = therefore.TimeLoop(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, theta)

therefore.MovieMaker(scalar_flux, L)

co = ['-k','-r','-b','-g','-m','-y', '-c', '-k','-r','-b','-g','-m','-y', '-c', '-k','-r','-b','-g','-m','-y', '-c','-k','-r','-b','-g','-m','-y', '-c']
print(co[7])

f=1
X = np.linspace(0, L, int(N_mesh*2+1))
x = np.linspace(0, L, int(N_mesh*2))
plt.figure(f)
plt.plot(x, np.sum(inital_angular_flux, axis=0), '-k')
for t in range(N_time):
    #ana = 2*analiticalSoultion(dt*(t+1), xsec, in_an_flux, source_mat, 1)
    #flatLinePlot(X, scalar_flux[:,t], co[t])
    plt.plot(x, scalar_flux[:,t])
    #plt.plot(5, ana, '^'+co[t])
    #flatLinePlot(X[N_mesh-10:N_mesh+11], scalar_flux[N_mesh-10:N_mesh+10,t], co[t])
plt.title('Infinte Med')
plt.xlabel('Distance')
plt.ylabel('Scalar Flux')
#plt.ylim([0,1.25*np.max(scalar_flux)]) #1.25*np.max(scalar_flux)
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
