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

def HasItConverged(scalar_flux_next, scalar_flux, tol=1e-6):
   error = max(abs((scalar_flux_next - scalar_flux) / scalar_flux_next))
   return(error < tol, error)
   
def RHSTransport(scalar_flux, scattering_xsec, source, N_mesh, dx):
    Q = np.zeros(N_mesh*2)
    for i in range(N_mesh):
        Q[2*i]   = scalar_flux[2*i]   * scattering_xsec[i]/2 + (source[i]/2)
        Q[2*i+1] = scalar_flux[2*i+1] * scattering_xsec[i]/2 + (source[i]/2)
    return(Q)
