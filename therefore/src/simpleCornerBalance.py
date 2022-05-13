import numpy as np
import numba as nb

#Simple Corner balence sweep
def SCBRun(angular_flux, Q, xsec, dx, mu, BCl, BCr, N_mesh):
    
    for angle in range(mu.size):
        if mu[angle] < 0: #goin back
            for i in range(N_mesh-1, -1, -1):
                #check bound
                if i == N_mesh-1:
                    psi_mh = BCr[angle] #
                else:
                    psi_mh = angular_flux[angle, 2*(i+1)]
                [angular_flux[angle, 2*i+1], angular_flux[angle, 2*i]] = SCBKernel_Linalg(Q[2*i+1], Q[2*i], psi_mh, xsec[i], dx[i], mu[angle])
        else: #goin forward
            for i in range(N_mesh):
                if i == 0:
                    psi_mh = BCl[angle] #
                else:
                    psi_mh = angular_flux[angle, 2*(i-1)+1]
                [angular_flux[angle, 2*i], angular_flux[angle, 2*i+1]] = SCBKernel_Linalg(Q[2*i], Q[2*i+1], psi_mh, xsec[i], dx[i], mu[angle])
                
    return(angular_flux)
    

#simple corner balance for a single cell
def SCBKernel(Q_entering, Q_exiting, psi_mh, xsec, dx, mu):
    
    mannaz = mu/2 + xsec*dx/2
    
    denominator = mannaz + (mu)**2/(4*mannaz)
    
    psi_entering = (dx/2*Q_entering - (mu*dx)/(4*mannaz)*Q_exiting + mu*psi_mh) / denominator
    
    psi_exiting = ((Q_exiting*dx)/2 + (mu*psi_entering)/2) / mannaz
    
    return(psi_entering, psi_exiting)


def SCBKernel_Linalg(Q_entering, Q_exiting, psi_mh, xsec, dx, mu):
    
    mannaz = mu/2 + xsec*dx/2
    
    A = np.array([[mannaz, -mu/2], [mu/2, mannaz]])
    
    b = np.array([[dx/2 * Q_exiting],[dx/2*Q_entering + mu*psi_mh]])
    
    [psi_exiting, psi_entering] = np.linalg.solve(A,b)
    
    return(psi_entering, psi_exiting)
