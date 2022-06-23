import numpy as np
from scipy.integrate import quad

#thanks Aaron!

def azurv1(x, c, t):
    '''Produce the anaylitical soultion for a vector
    returns a vector
    '''
    N_mesh = x.size
    
    scalar_flux = np.zeros(N_mesh)
    
    for i in range(N_mesh):
        scalar_flux[i] = phi(x[i], c, t)
    
    return(scalar_flux)
    
    
def azurv1_spav(x, c, t):
    '''Spatially averaged! x vector is location of cell bounds
    '''
    N_mesh = x.size-1
    dx = x[1] - x[0]
    
    scalar_flux = np.zeros(N_mesh)
    
    for i in range(N_mesh):
        scalar_flux[i] = quad(phi, x[i], x[i+1], args=(t, c), limit=1000)[0] / dx
    
    return(scalar_flux)

# Evaluate flux at position x and time t fir a scatterubg ratio c
def phi(x, c, t):
    
    if t == 0.0 or abs(x) >= t:
        return 0.0
        
    eta = x / t
    
    if c > 0:
        integral = quad(integrand,0,np.pi,args=(eta,t,c),limit=1000)[0]
    else:
        integral = 0
        
    return (np.exp(-t) / (2 * t)) *(1 + (c * t / (4 * np.pi)) * (1 - eta ** 2) * integral)

# Helper function for evaluating flux
def integrand(u,eta,t,c):
    i = complex(0,1)
    
    q   = (1+eta)/(1-eta)
    
    xi = (np.log(q)+i*u)/(eta+i*np.tan(u/2))
    return (1.0 / (np.cos(u / 2))) ** 2 * (xi ** 2 * np.exp(((c * t) / 2) * (1 - eta ** 2) * xi)).real
    
    

