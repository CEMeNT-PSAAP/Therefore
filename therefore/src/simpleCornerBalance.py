import numpy as np
import numba as nb

#Simple Corner balence sweep
#@nb.jit(nopython=True)
def SCBRun(angular_flux, Q, xsec, dx, mu, BCl, BCr, N_mesh):
    
    for angle in range(mu.size):
        if mu[angle] < 0: #goin back
            for i in range(N_mesh-1, -1, -1):
                #check bound
                if i == N_mesh-1:
                    psi_ph = BCr[angle] #
                else:
                    psi_ph = angular_flux[angle, 2*(i+1)]
                
                [angular_flux[angle, 2*i], angular_flux[angle, 2*i+1]]  = SCBKernel_Linalg_rtol(Q[2*i], Q[2*i+1], psi_ph, xsec[i], dx[i], mu[angle]) 
                
        else: #goin forward
            for i in range(N_mesh):
                
                if i == 0:
                    psi_mh = BCl[angle] #
                else:
                    psi_mh = angular_flux[angle, 2*(i-1)+1]
                
                [angular_flux[angle, 2*i], angular_flux[angle, 2*i+1]] = SCBKernel_Linalg_ltor(Q[2*i], Q[2*i+1], psi_mh, xsec[i], dx[i], mu[angle])
                
                
    return(angular_flux)
    
@nb.njit
def SCBKernel_Linalg_ltor(Q_r, Q_l, psi_mh, xsec, dx, mu):
    
    mannaz = mu/2 + xsec*dx/2
    
    A = np.array([[mannaz, mu/2],
                  [-mu/2, mannaz]])
    
    b = np.array([[dx/2*Q_l  + mu*psi_mh],
                  [dx/2*Q_r]])
    
    [psi_l, psi_r] = np.linalg.solve(A,b)
    
    return(psi_l, psi_r)

@nb.njit
def SCBKernel_Linalg_rtol(Q_r, Q_l, psi_ph, xsec, dx, mu):
    mannaz = xsec*dx/2 - mu/2 
    
    A = np.array([[mannaz, mu/2],
                  [-mu/2, mannaz]])
    
    b = np.array([[dx/2*Q_l],
                  [dx/2*Q_r - mu*psi_ph]])
    
    [psi_l, psi_r] = np.linalg.solve(A,b)
    
    return(psi_l, psi_r)
