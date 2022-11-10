from logging.handlers import QueueListener
import numpy as np
from .matrix import A_pos, A_neg, b_neg, b_pos, scatter_source
import therefore.src.utilities as utl
import numba as nb
np.set_printoptions(threshold=9999999)


def OCIMBSITimeStep(sim_perams, angular_flux_previous, source_mesh, xsec_mesh, xsec_scatter_mesh, dx_mesh, angles, weights):

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

    angular_flux = np.zeros([N_angles, N_ans], data_type)
    angular_flux_mid = np.zeros([N_angles, N_ans], data_type)

    angular_flux_last = np.zeros([N_angles, N_ans], data_type)
    angular_flux_mid_last = np.zeros([N_angles, N_ans], data_type)

    scalar_flux = np.zeros(N_ans, data_type)
    scalar_flux_mid = np.zeros(N_ans, data_type)

    scalar_flux_last = np.zeros(N_ans, data_type)
    scalar_flux_next = np.zeros(N_ans, data_type)

    scalar_flux_mid_next = np.zeros(N_ans, data_type)

    while source_converged == False:

        if printer:
            print()

            #print()
            #print(scalar_flux_next)
            #print()
            #print(scalar_flux_mid_next)

        BCl = utl.BoundaryCondition(sim_perams['boundary_condition_left'],   0, N_mesh, angular_flux=angular_flux, incident_flux_mag=sim_perams['left_in_mag'],  angle=sim_perams['left_in_angle'],  angles=angles)
        BCr = utl.BoundaryCondition(sim_perams['boundary_condition_right'], -1, N_mesh, angular_flux=angular_flux, incident_flux_mag=sim_perams['right_in_mag'], angle=sim_perams['right_in_angle'], angles=angles)

        x = 0
        [angular_flux, angular_flux_mid] = Itteration(angular_flux_previous, angular_flux_last, angular_flux_mid_last, scalar_flux, scalar_flux_mid_next, source_mesh, xsec_mesh, xsec_scatter_mesh, dx_mesh, dt, velocity, angles, BCl, BCr)

        #calculate current
        current = utl.Current(angular_flux, weights, angles)
        
        #calculate scalar flux for next itteration
        scalar_flux_next = utl.ScalarFlux(angular_flux, weights)
        scalar_flux_mid_next = utl.ScalarFlux(angular_flux_mid, weights)

        #print()
        #print(scalar_flux_next)
        #print()
        #print(scalar_flux_mid_next)

        if source_counter > 2:
            #check for convergence
            error_eos = np.linalg.norm(scalar_flux_mid_next - scalar_flux_mid, ord=2)
            error_mos = np.linalg.norm(scalar_flux_next - scalar_flux, ord=2)

            #print(error_eos)
            #print(error_mos)
            #print()
            #print()

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


def Itteration(angular_flux_previous, angular_flux_last, angular_flux_midstep_last, scalar_flux, scalar_flux_halfNext, Q, xsec, xsec_scatter, dx, dt, v, mu, BCl, BCr):
    N_angle = mu.size
    N_mesh = dx.size

    sizer: int = 4

    half = int(mu.size/2)

    angular_flux_next = np.zeros_like(angular_flux_previous)
    angular_flux_mid_next = np.zeros_like(angular_flux_previous)

    for angle in range(N_angle):
        for i in range(N_mesh):
            i_l: int = int(2*i)
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

            A = np.zeros((sizer,sizer))
            b = np.zeros((sizer,1))

            if mu[angle] < 0:
                A = A_neg(dx[i], v, dt, mu[angle], xsec[i])
                b = b_neg(dx[i], v, dt, mu[angle], Ql, Qr, Ql, Qr, psi_halfLast_L, psi_halfLast_R, psi_rightBound, psi_halfNext_rightBound, xsec_scatter[i], scalar_flux[i_l], scalar_flux[i_r], scalar_flux_halfNext[i_l], scalar_flux_halfNext[i_r])
                #b_neg(dx, v, dt, mu, Ql, Qr, Q_halfNext_L, Q_halfNext_R, psi_halfLast_L, psi_halfLast_R, psi_rightBound, psi_halfNext_rightBound, xsec_scatter, phi_L, phi_R, phi_halfNext_L, phi_halfNext_R):

            elif mu[angle] > 0:
                A = A_pos(dx[i], v, dt, mu[angle], xsec[i])
                b = b_pos(dx[i], v, dt, mu[angle], Ql, Qr, Ql, Qr, psi_halfLast_L, psi_halfLast_R, psi_leftBound, psi_halfNext_leftBound, xsec_scatter[i], scalar_flux[i_l], scalar_flux[i_r], scalar_flux_halfNext[i_l], scalar_flux_halfNext[i_r])
            
            psi_raw = np.linalg.solve(A, b)

            angular_flux_next[angle, i_l] = psi_raw[0]
            angular_flux_next[angle, i_r] = psi_raw[1]
            angular_flux_mid_next[angle, i_l] = psi_raw[2]
            angular_flux_mid_next[angle, i_r] = psi_raw[3]

    return(angular_flux_next, angular_flux_mid_next)




#def cell():
    


