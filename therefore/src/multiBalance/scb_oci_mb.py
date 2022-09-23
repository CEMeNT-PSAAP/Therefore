import numpy as np
from timeTest import N_angle
from .matrix import A_pos, A_neg, c_neg, c_pos, scatter_source

# last refers to iteration
# prev refers to previous time step
def OCIMBRun(angular_flux_last, angular_flux_midstep_last, angular_flux_midstep_previous, source, xsec, xsec_scatter, dx, dt, v, mu, weight, BCl, BCr):
    
    N_mesh = int(dx.size)
    angular_flux = np.zeros_like(angular_flux_last)
    angular_flux_midstep = np.zeros_like(angular_flux_last)
    half = int(mu.size/2)

    
    sizer = int(mu.size*4)
    N_angle = mu.size

    A = np.zeros([sizer,sizer])
    c = np.zeros([sizer,1])

    for i in range(N_mesh):
        
        i_l = int(2*i)
        i_r = int(2*i+1)

        Q = source[:,i_l:i_r]

        psi_halfLast_L = angular_flux_midstep_previous[:, i_l] # known all angles from the last time step in the left cell
        psi_halfLast_R = angular_flux_midstep_previous[:, i_r] # known

        # getting really sloppy with the indiceis
        for m in range(N_angle):
            # boundary conditions
            if i == 0:  #left
                psi_rightBound = angular_flux_last[m, i_r+1] # iterating on (unknown)
                psi_leftBound =  BCl[m] # known

                psi_halfNext_rightBound = angular_flux_midstep_last[m, i_r+1] # iterating on (unknown)
                psi_halfNext_leftBound  = BCl[m] # known

            elif i == N_mesh-1: #right
                psi_rightBound = BCr[m] # known
                psi_leftBound =  angular_flux_last[m, i_l-1] # iterating on (unknown)

                psi_halfNext_rightBound = BCr[m] # known
                psi_halfNext_leftBound  = angular_flux_midstep_last[m, i_l-1] # iterating on (unknown)

            else: #middles
                psi_rightBound = angular_flux_last[m, i_r+1] # iterating on (unknown)
                psi_leftBound =  angular_flux_last[m, i_l-1] # iterating on (unknown)

                psi_halfNext_rightBound = angular_flux_midstep_last[m, i_r+1] # iterating on (unknown)
                psi_halfNext_leftBound  = angular_flux_midstep_last[m, i_l-1] # iterating on (unknown)


            if mu[m] < 0:
                A_small = A_neg(dx, v, dt, mu[m], xsec)
                c_small = c_neg(dx, v, dt, mu[m], Q[m,i_l], Q[m,i_r], Q[m,i_l], Q[m,i_r], psi_halfLast_L[m], psi_halfLast_R[m], psi_rightBound, psi_halfNext_rightBound)

            elif mu[m] > 0:
                A_small = A_pos(dx, v, dt, mu[m], xsec)
                c_small = c_pos(dx, v, dt, mu[m], Q[m,i_l], Q[m,i_r], Q[m,i_l], Q[m,i_r], psi_halfLast_L[m], psi_halfLast_R[m], psi_leftBound, psi_halfNext_leftBound)

            S_small = scatter_source(dx, xsec_scatter, N_angle, weight[i])

            A[m*4:(m+1)*4, m*4:(m+1)*4] = A_small - S_small
            c[m*4:(m+1)*4] = c_small
            
        angular_flux_raw = np.linalg.solve(A, c)

        for m in range(N_angle):
                angular_flux[m,i_l]         = angular_flux_raw[4*m]
                angular_flux[m,i_r]         = angular_flux_raw[4*m+1]
                
                angular_flux_midstep[m,i_l] = angular_flux_raw[4*m+2]
                angular_flux_midstep[m,i_r] = angular_flux_raw[4*m+3]

    
    return(angular_flux, angular_flux_midstep)
