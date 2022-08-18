import numpy as np


#utility functions
def ScalarFlux(angular_flux, weights):
    scalar_flux = np.zeros(angular_flux.shape[1])
    for i in range(angular_flux.shape[0]):
        scalar_flux += weights[i] * angular_flux[i,:]
        
    return(scalar_flux)

def Current(angular_flux, angles, weights):
    current = np.zeros(angular_flux.shape[1])
    for i in range(angular_flux.shape[0]):
        current += weights[i] * angles[i] * angular_flux[i,:]
    return(current)

def HasItConverged(a, b, tol=1e-16):
   close = np.allclose(a, b, atol=tol)
   return(close)
   
def RHSTransport(scalar_flux, scattering_xsec, source, N_mesh, dx):
    N = 6
    #N_angle = source.shape[0]
    #print(source.shape)
    #print(scalar_flux.shape)
    Q = np.zeros([N])
    print(N)
    for i in range(N):
        #for j in range(N_angle):
        Q[i] = scalar_flux[i] * scattering_xsec[int(i/2)]/2 + (source[int(i/2)])
    return(Q)
    
