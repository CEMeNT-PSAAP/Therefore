import numpy as np
import numba as nb

#@nb.jit(nopython=True, parallel=True)
def OCIRun(angular_flux, source, xsec, xsec_scatter, dx, mu, weight, BCl, BCr):
    
    n_mesh = int(dx.size)
    angular_flux_next = np.zeros_like(angular_flux)
    
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



#@nb.jit(nopython=True)
def SCB_OneCellInv_Cell(angular_flux, source, xsec, xsec_scatter, dx, mu, weight):
    
    n_angle = mu.size
    
    A = np.zeros((n_angle*2, n_angle*2), dtype=np.float64)
    b = np.zeros(n_angle*2, dtype=np.float64)
    
    const = -(xsec_scatter*dx) / 2
    
    alpha = xsec*dx/2
    #build matrix
    #lefts [k,i], [row, col] so we are filling in col by col
    for k in range (0, n_angle):
        beta  = mu[k]/2
        
        #lefts
        A[2*k,k*2]   = alpha + beta
        A[2*k,k*2+1] = beta
        
        #rights
        A[2*k+1,k*2]   = -beta
        A[2*k+1,k*2+1] = -beta + alpha
        
        #scalar flux
        for i in range (0, n_angle):
            A[2*k,2*i]     += const*weight[i]
            A[2*k+1,2*i+1] += const*weight[i]
            
    for i in range (n_angle):
        b[2*i]   = (source*dx)/2 + mu[i] * angular_flux[i, 0]
        b[2*i+1] = (source*dx)/2 - mu[i] * angular_flux[i, 1]

    #print(A)
    #print()
    #print(b)    
    #print()
    #print('Spec rad: {0}'.format(max(abs(np.linalg.eigvals(A)))))
    
    
    next_angflux = np.linalg.solve(A, b).reshape((-1, 2))
    
    #print()
    #print(next_angflux)
    
    return(next_angflux)



def neg_flux_fixup(next_angflux):
    for i in range(next_angflux.shape[0]):
        for k in range(next_angflux.shape[1]):
            if next_angflux[i,k] < 0:
                next_angflux[i,k] = 0
    
if __name__ == '__main__':
    
    
    angular_flux  = np.array([[1,1],[1,1]])
    xsec = 1
    xsec_scatter = 0
    dx = 2
    mu = np.array([1,1])
    weights = np.array([1,1])
    N_mesh = 0
    source = 1
    
    SCB_OneCellInv_Cell(angular_flux, source, xsec, xsec_scatter, dx, mu, weights)
