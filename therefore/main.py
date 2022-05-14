"""
Therefore
Method of implementation for slab wall problem
Author: Jackson P. Morgan
breif: An implementaiton of simplified corner balance for assignemtn 2 in NSE 653
prof: Todd Palmer
date: May 9th 2022
"""

import numpy as np
import matplotlib.pyplot as plt
import numba as nb
import therefore.src as src

"""
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
        

def BoundaryCondition(BC, i, N_mesh, angular_flux=None, incident_flux_mag=None, angles=None, angle=None):

    if BC == 'reflecting':
        N_angles: int = angular_flux.shape[0]
        half = int(N_angles/2)
        psi_required = np.zeros(N_angles)
        
        if i == 0: #left edge
            psi_required[half:] = angular_flux[:half, 0]
            
        else:
            psi_required[:half] = angular_flux[half:, -1]

        
    elif BC == 'vacuum':
        psi_required = np.zeros(angular_flux.shape[0])
    
    elif BC == 'incident_iso':
        psi_required = BC_isotropic(incident_flux_mag, angles)
        
    elif BC == 'incident_ani':
        psi_required = BC_ani(incident_flux_mag, angle, angles)
        
    else:
        print()
        print('>>>Error: No Boundary Condition Supplied<<<')
        print()
    
    return(psi_required)

def BC_isotropic(incident_flux_mag, angles):
    BC = (incident_flux_mag/angles.size)*np.ones(angles.size)
    return(BC)

def BC_ani(incident_flux_mag, angle, angles):
    angle_id = find_nearest(angles, angle)
    BC = np.zeros(angles.size)
    BC[angle_id] = incident_flux_mag
    return(BC)
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return (idx)

def RHS_transport(scalar_flux, scattering_xsec, source, N_mesh, dx):
    Q = np.zeros(N_mesh*2)
    for i in range(N_mesh):
        Q[2*i]   = scalar_flux[2*i]   * scattering_xsec[i]/2 + (source[i]/2)
        Q[2*i+1] = scalar_flux[2*i+1] * scattering_xsec[i]/2 + (source[i]/2)
    return(Q)
"""

def SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh):
    
        #print title contents
    #with open("title_print.txt", "r", encoding="utf-8") as file:
    #    for line in file:
    #        print(line.strip())
    
    #Inputs
    order_gauss_quad = sim_perams['N_angles']
    N_angles = order_gauss_quad
    
    data_type = sim_perams['data_type']
    
    
    
    boundary_condition_left =  sim_perams['boundary_condition_left'] #incident_iso
    boundary_condition_right = sim_perams['boundary_condition_right'] #'reflecting'
    left_in_mag = sim_perams['left_in_mag']
    right_in_mag = sim_perams['right_in_mag']
    left_in_angle = sim_perams['left_in_angle']
    right_in_angle = sim_perams['right_in_angle']
    
    #mesh building
    L = sim_perams['L']
    N_mesh = sim_perams['N_mesh']
    
    
    
    
    """
    #mesh builder for Reed's Problem
    region_id = np.array([1,2,3,4,5], int)
    region_widths = np.array([2,2,1,1,2], int)
    region_bounds = np.array([2,4,5,6,8], float)
    sigma_s = np.array([.9, .9, 0, 0, 0], data_type)
    sigma_t = np.array([1, 1, 0, 5, 50], data_type)
    Source = np.array([0, 1, 0, 0, 50], data_type)
    dx = .1/sigma_t
    dx[2] = .25
    N_region = np.array(region_widths/dx, int)
    
    N_mesh: int = sum(N_region)
    
    xsec_mesh = np.empty(N_mesh, data_type)
    xsec_scatter_mesh = np.empty(N_mesh, data_type)
    dx_mesh = np.empty(N_mesh, data_type)
    source_mesh = np.empty(N_mesh, data_type)
    region_id_mesh = np.empty(N_mesh, data_type)
    
    for i in range(region_widths.size):
        LB = sum(N_region[:i])
        RB = sum(N_region[:i+1])
        xsec_mesh[LB:RB] = sigma_t[i]
        xsec_scatter_mesh[LB:RB] = sigma_s[i]
        dx_mesh[LB:RB] = dx[i]
        source_mesh[LB:RB] = Source[i]
        region_id_mesh[LB:RB] = region_id[i]
        
        
    boundary_condition_left =  'vacuum' #incident_iso
    boundary_condition_right = 'reflecting' #'reflecting'
    
    # END OF REEDS PROBLEM
    """
    [angles_gq, weights_gq] = np.polynomial.legendre.leggauss(order_gauss_quad)

    angular_flux = np.zeros([order_gauss_quad, int(N_mesh*2)], data_type)
    
    scalar_flux  = np.zeros(int(N_mesh*2), data_type)
    scalar_flux_next  = np.zeros(int(N_mesh*2), data_type)


    # TODO: Source itterations
    source_converged = False
    source_counter = 0
    
    while source_converged == False:
        
        print('Next Itteration: {0}'.format(source_counter),end='\r')
        
        BCl = src.BoundaryCondition(boundary_condition_left,   0, N_mesh, angular_flux=angular_flux, incident_flux_mag=left_in_mag, angle=left_in_angle, angles=angles_gq)
        BCr = src.BoundaryCondition(boundary_condition_right, -1, N_mesh, angular_flux=angular_flux, incident_flux_mag=right_in_mag, angle=right_in_angle, angles=angles_gq)
        
        Q = src.RHSTransport(scalar_flux, xsec_scatter_mesh, source_mesh, N_mesh, dx_mesh)
        
        
        # TODO: simple corner balance
        angular_flux = src.SCBRun(angular_flux, Q, xsec_mesh, dx_mesh, angles_gq, BCl, BCr, N_mesh)
        
        # TODO: calculate current
        current = src.Current(angular_flux, weights_gq, angles_gq)
        
        # TODO: calculate scalar flux for next itteration
        scalar_flux_next = src.ScalarFlux(angular_flux, weights_gq)
        
        
        # TODO: Check for convergence
        source_converged, error = src.HasItConverged(scalar_flux_next, scalar_flux)
        
        if source_counter > 10000:
            print('Error source not converged after 1000 itterations')
            print()
            source_converged = True
        
        scalar_flux = scalar_flux_next
        source_counter += 1
    
    return(scalar_flux, current)
    
    
if __name__ == '__main__':
    main()
