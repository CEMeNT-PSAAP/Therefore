import numpy as np


#utility functions
def ScalarFlux(angular_flux, weights):
    scalar_flux = np.zeros(angular_flux.shape[1])
    for i in range(angular_flux.shape[0]):
        #print('i in scalar flux comp {0}'.format(i))
        scalar_flux += weights[i] * angular_flux[i,:]
        
    return(scalar_flux)

def Current(angular_flux, angles, weights):
    current = np.zeros(angular_flux.shape[1])
    for i in range(angular_flux.shape[0]):
        current += weights[i] * angles[i] * angular_flux[i,:]
    return(current)

def HasItConverged(scalar_flux_next, scalar_flux, tol=1e-8):
   error = max(abs((scalar_flux_next - scalar_flux) / scalar_flux_next))
   return(error < tol, error)
