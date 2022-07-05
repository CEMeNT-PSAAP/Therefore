import numpy as np
import numba as nb
import scipy.linalg as sci

#@nb.jit(nopython=True, parallel=True)
def OCIRun(angular_flux, source, xsec, xsec_scatter, dx, mu, weight, BCl, BCr):
    
    n_mesh = int(dx.size)
    angular_flux_next = np.zeros_like(angular_flux)
    half = int(mu.size/2)
    #print(n_mesh)
    for i in nb.prange(n_mesh):
        
        bound_ang_flux = np.zeros((mu.size, 2), dtype=np.float64)
        Q_flux = np.zeros((mu.size, 2), dtype=np.float64)
        Q_flux = source[:,2*i:2*i+2]
        
        if i == n_mesh-1:   #RHS BC
            bound_ang_flux[:,0] = angular_flux[:,-3]
            bound_ang_flux[:,1] = BCr
        elif i == 0:        #LHS BC
            bound_ang_flux[:,0] = BCl
            bound_ang_flux[:,1] = angular_flux[:,2]
        else:               #interior cell
            bound_ang_flux[:,0] = angular_flux[:,i*2-1]
            bound_ang_flux[:,1] = angular_flux[:,i*2+2]
            
        #bound_ang_flux[0,0] = 20
        #bound_ang_flux[1,1] = 20
            
        angular_flux_cell = SCB_OneCellInv_Cell(bound_ang_flux, Q_flux, xsec[i], xsec_scatter[i], dx[i], mu, weight)
        
        angular_flux_next[:,2*i] = angular_flux_cell[:,0]
        angular_flux_next[:,2*i+1] = angular_flux_cell[:,1]
    
    return(angular_flux_next)



#@nb.jit(nopython=True)
def SCB_OneCellInv_Cell(angular_flux, source, xsec, xsec_scatter, dx, mu, weight):
    
    n_angle = mu.size
    half = int(n_angle/2)
    
    A = np.zeros((n_angle*2, n_angle*2), dtype=np.float64)
    b = np.zeros(n_angle*2, dtype=np.float64)
    
    const = (xsec_scatter*dx)/4
    
    alpha = xsec*dx/2
    
    #negative ordinants
    #'''
    for i in range(half): 
        A[2*i, 2*i] = -mu[i]/2
        A[2*i, 2*i+1] = -mu[i]/2 + alpha
        A[2*i+1, 2*i] = -mu[i]/2 + alpha
        A[2*i+1, 2*i+1] = mu[i]/2
        
        b[2*i] =   (source[i,0]*dx)/2 - (mu[i] * angular_flux[i, 1])
        b[2*i+1] = (source[i,1]*dx)/2
    
    #positive ordinants
    for i in range(half, n_angle, 1): 
        A[2*i, 2*i] = mu[i]/2 + alpha
        A[2*i, 2*i+1] = mu[i]/2
        A[2*i+1, 2*i] = -mu[i]/2
        A[2*i+1, 2*i+1] = mu[i]/2 + alpha
    	
        b[2*i] =   (source[i,0]*dx)/2 + (mu[i] * angular_flux[i, 0])
        b[2*i+1] = (source[i,1]*dx)/2
        #'''
    #scalar flux
    for k in range (0, n_angle):
        for i in range (0, n_angle):
            A[2*k,2*i]     -= const*weight[i]
            A[2*k+1,2*i+1] -= const*weight[i]
    
    print(A)
    print()
    print()
    #print(b)
    
    next_angflux = np.linalg.solve(A, b).reshape((-1, 2))
    
    return(next_angflux)



def neg_flux_fixup(next_angflux):
    for i in range(next_angflux.shape[0]):
        for k in range(next_angflux.shape[1]):
            if next_angflux[i,k] < 0:
                next_angflux[i,k] = 0
    
if __name__ == '__main__':
    angular_flux  = np.array([[1,0],[0,0]])
    xsec = 1
    xsec_scatter = 1
    dx = 2
    mu = np.array([-.235,.235])
    weights = np.array([1,1])
    N_mesh = 0
    source = np.array([[0,0],[0,0]])
    
    a_out = SCB_OneCellInv_Cell(angular_flux, source, xsec, xsec_scatter, dx, mu, weights)
    
    print(a_out)
    
    
