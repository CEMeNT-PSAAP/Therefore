from logging.handlers import QueueListener
import numpy as np
from .matrix import A_pos, A_neg, b_neg, b_pos, scatter_source
import therefore.src.utilities as utl
import numba as nb

import cupyx.scipy.sparse.linalg as gpuLinalg
import scipy.sparse.linalg as cpuLinalg
import cupyx.scipy.sparse as spMat
import cupy as cu
from scipy.sparse import csr_matrix, lil_matrix
import betterspy

np.set_printoptions(threshold=9999999)


def SIMBTimeStepBig(sim_perams, angular_flux_previous, angular_flux_mid_previous, source_mesh, xsec_mesh, xsec_scatter_mesh, dx_mesh, angles, weights):

    velocity = sim_perams['velocity']
    dt = sim_perams['dt']
    N_time = sim_perams['N_time']
    N_angles = sim_perams['N_angles']
    N_mesh = sim_perams['N_mesh']
    data_type = sim_perams['data_type']
    max_it = sim_perams['max loops']
    tol = sim_perams['tolerance']
    printer = sim_perams['print']

    source_converged: bool = False
    source_counter: int = 0
    no_convergence: bool = False
    spec_rad = 0
    
    N_ans = int(2*N_mesh)

    #source_mesh /= 2

    angular_flux = np.zeros([N_angles, N_ans], data_type)
    angular_flux_mid = np.zeros([N_angles, N_ans], data_type)

    angular_flux_last = np.zeros([N_angles, N_ans], data_type)
    angular_flux_mid_last = np.zeros([N_angles, N_ans], data_type)

    #angular_flux_last = angular_flux_previous ##
    #angular_flux_mid_last = angular_flux_mid_previous # #

    scalar_flux = np.zeros(N_ans, data_type)
    scalar_flux_mid = np.zeros(N_ans, data_type)

    scalar_flux_last = np.zeros(N_ans, data_type)
    scalar_flux_next = np.zeros(N_ans, data_type)

    scalar_flux_mid_next = np.zeros(N_ans, data_type)

    A = buildHer(xsec_mesh, dx_mesh, dt, velocity, angles)


    while source_converged == False:

        BCl = utl.BoundaryCondition(sim_perams['boundary_condition_left'],   0, N_mesh, angular_flux=angular_flux, incident_flux_mag=sim_perams['left_in_mag'],  angle=sim_perams['left_in_angle'],  angles=angles)
        BCr = utl.BoundaryCondition(sim_perams['boundary_condition_right'], -1, N_mesh, angular_flux=angular_flux, incident_flux_mag=sim_perams['right_in_mag'], angle=sim_perams['right_in_angle'], angles=angles)

        b = buildHim(angular_flux_previous, angular_flux_last, angular_flux_mid_last, scalar_flux, scalar_flux_mid_next, source_mesh, xsec_scatter_mesh, dx_mesh, dt, velocity, angles, BCl, BCr)
        angular_flux_raw, info = cpuLinalg.gmres(A, b)

        print(np.transpose(angular_flux_raw))

        [angular_flux, angular_flux_mid] = resort(angular_flux_raw, angles, N_mesh)

        #calculate current
        current = utl.Current(angular_flux, weights, angles)
        
        #calculate scalar flux for next itteration
        scalar_flux_next = utl.ScalarFlux(angular_flux, weights)
        scalar_flux_mid_next = utl.ScalarFlux(angular_flux_mid, weights)

        if source_counter > 2:
            #check for convergence
            error_eos = np.linalg.norm(scalar_flux_mid_next - scalar_flux_mid, ord=2)
            error_mos = np.linalg.norm(scalar_flux_next - scalar_flux, ord=2)

            if error_eos < tol and error_mos < tol:
                source_converged = True

            spec_rad = np.linalg.norm(scalar_flux_next - scalar_flux, ord=2) / np.linalg.norm((scalar_flux - scalar_flux_last), ord=2)

        if source_counter > max_it:
            print('Error source not converged after max iterations')
            print()
            source_converged = True
            no_convergence = True

        angular_flux_last = angular_flux
        angular_flux_mid_last = angular_flux_mid
        
        scalar_flux_last = scalar_flux
        scalar_flux = scalar_flux_next

        scalar_flux_mid = scalar_flux_mid_next

        source_counter += 1

    return(angular_flux, angular_flux_mid, current, spec_rad, source_counter, source_converged)



def buildHer(xsec, dx, dt, v, mu):
    N_angle = mu.size
    N_mesh = dx.size

    sizer: int = 4 * N_mesh

    A_uge = lil_matrix((4*N_angle*N_mesh, 4*N_angle*N_mesh))

    for j in range(N_angle):
        A_quadrature = lil_matrix((sizer, sizer))


        for i in range(N_mesh):
            # forming the lower triangular mats 
            # cols of our A_quadrature matrix are filled with the same sub mats
            # forming lower triangular for positive sweeps and reflected for negaitive sweeps

            start = 0
            stop = i+1
            step = 1

            if   mu[j] < 0: # negaitive
                A = A_neg(dx[i], v, dt, mu[j], xsec[i])
            elif mu[j] > 0: # positive
                A = A_pos(dx[i], v, dt, mu[j], xsec[i])
            
            for k in range(start, stop, step):
                row_left:  int = 4*i
                row_right: int = 4*(i+1)
                col_top:   int = 4*k
                col_bot:   int = 4*(k+1)
                A_quadrature[row_left:row_right, col_top:col_bot] = A

        Ba = 4*N_mesh*j
        Bb = 4*N_mesh*(j+1)

        A_uge[Ba:Bb, Ba:Bb] = A_quadrature

    return(A_uge)


def buildHim(angular_flux_previous, angular_flux_last, angular_flux_midstep_last, scalar_flux, scalar_flux_halfNext, Q, xsec_scatter, dx, dt, v, mu, BCl, BCr):
    N_angle = mu.size
    N_mesh = dx.size

    sizer: int = 4 * N_angle

    b_big = np.zeros((4*N_angle*N_mesh))

    for angle in nb.prange(N_angle):
        for i in range(N_mesh):
            i_l: int = int(2*i) # left index of convencnence
            i_r: int = int(2*i+1)

            psi_halfLast_L = angular_flux_previous[angle, i_l]
            psi_halfLast_R = angular_flux_previous[angle, i_r]

            Ql = Q[angle, i_l]
            Qr = Q[angle, i_r]

            if i == 0:  #left
                psi_rightBound = angular_flux_last[angle, i_r+1] # iterating on (unknown)
                psi_leftBound =  BCl[angle] # known

                psi_halfNext_rightBound = angular_flux_midstep_last[angle, i_r+1] # iterating on (unknown)
                psi_halfNext_leftBound  = BCl[angle] # known

            elif i == N_mesh-1: #right
                psi_rightBound = BCr[angle] # known
                psi_leftBound =  angular_flux_last[angle, i_l-1] # iterating on (unknown)

                psi_halfNext_rightBound = BCr[angle] # known
                psi_halfNext_leftBound  = angular_flux_midstep_last[angle, i_l-1] # iterating on (unknown)

            else: #middles
                psi_rightBound = angular_flux_last[angle, i_r+1] # iterating on (unknown)
                psi_leftBound =  angular_flux_last[angle, i_l-1] # iterating on (unknown)

                psi_halfNext_rightBound = angular_flux_midstep_last[angle, i_r+1] # iterating on (unknown)
                psi_halfNext_leftBound  = angular_flux_midstep_last[angle, i_l-1] # iterating on (unknown)

            if mu[angle] < 0:
                b = b_neg(dx[i], v, dt, mu[angle], Ql, Qr, Ql, Qr, psi_halfLast_L, psi_halfLast_R, psi_rightBound, psi_halfNext_rightBound, xsec_scatter[i], scalar_flux[i_l], scalar_flux[i_r], scalar_flux_halfNext[i_l], scalar_flux_halfNext[i_r])

            elif mu[angle] > 0:
                b = b_pos(dx[i], v, dt, mu[angle], Ql, Qr, Ql, Qr, psi_halfLast_L, psi_halfLast_R, psi_leftBound, psi_halfNext_leftBound, xsec_scatter[i], scalar_flux[i_l], scalar_flux[i_r], scalar_flux_halfNext[i_l], scalar_flux_halfNext[i_r])
            
            Bup = angle*N_mesh + 4*i
            Bbo = angle*N_mesh + 4*(i+1)

            b_big[Bup:Bbo] = b[:,0]

    return(b_big)


def resort(angular_flux_raw, angles, N_mesh):
    N_angles = angles.size

    af = np.zeros((N_angles, N_mesh*2))
    af_mid = np.zeros((N_angles, N_mesh*2))

    for i in range(N_angles):
        #print()
        #print('angle {0}'.format(i))
        start_LU = int(i*N_mesh*4)
        #print('first element in that angle {0}'.format(start_LU))

        if angles[i] > 0:
            #print('Angle Positive!')
            for j in range(N_mesh):

                x = start_LU + j*4
                print('angle {0}'.format(i))
                print('first element in that angle {0}'.format(start_LU))
                print('for loop {0}, raw index {1}'.format(j,x))

                af[i,2*j]           = angular_flux_raw[x]
                af[i,2*j+1]         = angular_flux_raw[x+1]
                
                af_mid[i,2*j]   = angular_flux_raw[x+2]
                af_mid[i,2*j+1] = angular_flux_raw[x+3]


        elif angles[i] < 0:
            #print('Angle Negative!')
            for j in range(N_mesh):
                x = start_LU + (N_mesh - j)*4

                #print('for loop {0}, raw index {1}'.format(j,x))

                af[i,2*j]           = angular_flux_raw[x]
                af[i,2*j+1]         = angular_flux_raw[x+1]
                
                af_mid[i,2*j]   = angular_flux_raw[x+2]
                af_mid[i,2*j+1] = angular_flux_raw[x+3]

    return(af, af_mid)


if __name__ == '__main__':
    

    N_mesh = int(10)

    xsec = .5*np.ones(N_mesh)
    xsec_scatter = 2*np.ones(N_mesh)
    dx = .75*np.ones(N_mesh)
    dt = .5
    v = 1
    weight = np.array([1,1])
    mu = np.array([-1, -.5, .5, 1])

    A = buildHer(xsec, dx, dt, v, mu)

    betterspy.show(A)
    betterspy.write_png("out.png", A)


#def cell():