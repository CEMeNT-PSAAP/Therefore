import numpy as np
import numba as nb
import scipy.linalg as sci

@nb.jit(nopython=True, parallel=True)
def OCIRun(angular_flux, source, xsec, xsec_scatter, dx, mu, weight, BCl, BCr):
    
    n_mesh = int(dx.size)
    angular_flux_next = np.zeros_like(angular_flux)
    half = int(mu.size/2)
    
    for i in nb.prange(n_mesh):
        
        bound_ang_flux = np.zeros((mu.size, 2), dtype=np.float64)
        
        if i == n_mesh-1:   #RHS BC
            bound_ang_flux[:,0] = angular_flux[:,-3]
            bound_ang_flux[:,1] = BCr
        elif i == 0:        #LHS BC
            bound_ang_flux[:,0] = BCl
            bound_ang_flux[:,1] = angular_flux[:,2]
        else:               #interior cell
            bound_ang_flux[:,0] = angular_flux[:,i*2-1]
            bound_ang_flux[:,1] = angular_flux[:,i*2+2]
            
        angular_flux_cell = SCB_OneCellInv_Cell(bound_ang_flux, source[i], xsec[i], xsec_scatter[i], dx[i], mu, weight)
        
        angular_flux_next[:,2*i] = angular_flux_cell[:,0]
        angular_flux_next[:,2*i+1] = angular_flux_cell[:,1]

    
    return(angular_flux_next)



@nb.jit(nopython=True)
def SCB_OneCellInv_Cell(angular_flux, source, xsec, xsec_scatter, dx, mu, weight):
    
    n_angle = mu.size
    
    A = np.zeros((n_angle*2, n_angle*2), dtype=np.float64)
    b = np.zeros(n_angle*2, dtype=np.float64)
    
    const = -(xsec_scatter*dx)/4
    
    alpha = xsec*dx/2
    
    #build matrix
    #lefts [k,i], [row, col] so we are filling in col by col
    A[0,0] = -mu[0]/2
    A[0,1] = -mu[0]/2 + alpha
    
    A[1,0] = -mu[0]/2 + alpha
    A[1,1] =  mu[0]/2
    
    A[2,2] =  mu[1]/2 + alpha
    A[2,3] =  mu[1]/2
    
    A[3,2] = -mu[1]/2
    A[3,3] =  mu[1]/2 + alpha
    
    
    for k in range (0, n_angle):
        for i in range (0, n_angle):
            A[2*k,2*i]     += const*weight[i]
            A[2*k+1,2*i+1] += const*weight[i]
            
    b[0] = (source*dx)/2 - (mu[0] * angular_flux[0, 1])
    b[1] = (source*dx)/2
    b[2] = (source*dx)/2 + (mu[1] * angular_flux[1, 0])
    b[3] = (source*dx)/2
    
    
    
    next_angflux = np.linalg.solve(A, b).reshape((-1, 2))
    #next_angflux = next_angflux.reshape((-1, 2))
    #next_angflux = np.zeros([2,2])
    
    #next_angflux[0,0] = next_linalg_exp[0]
    #next_angflux[0,1] = next_linalg_exp[1]
    #next_angflux[1,0] = next_linalg_exp[2]
    #next_angflux[1,1] = next_linalg_exp[3]
    
    
    #print()
    #print(next_angflux)
    
    return(next_angflux)



def neg_flux_fixup(next_angflux):
    for i in range(next_angflux.shape[0]):
        for k in range(next_angflux.shape[1]):
            if next_angflux[i,k] < 0:
                next_angflux[i,k] = 0
    
if __name__ == '__main__':
    
    
    angular_flux  = np.array([[1,0],[0,0]])
    xsec = 1
    xsec_scatter = 0
    dx = 2
    mu = np.array([-.235,.235])
    weights = np.array([1,1])
    N_mesh = 0
    source = 0
    
    a_out = SCB_OneCellInv_Cell(angular_flux, source, xsec, xsec_scatter, dx, mu, weights)
    
    print(a_out)
    
    
