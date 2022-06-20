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
import os
import sys
import scipy as sci
import math

def SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, ang_flux_last=None, ret_angflux=False):
    ''' Return converged scalar flux and current
    
    Implementation of simple corner balance source itterations
    '''
    
    #if running from command line print title contents
    #with open(os.path.join(sys.path[0], "title_print.txt"), "r", encoding="utf-8") as file:
    #    for line in file:
    #        print(line.strip())
    
    #Inputs: sets inputs from sim_perams tuple
    order_gauss_quad = sim_perams['N_angles']
    N_angles = order_gauss_quad
    
    data_type = sim_perams['data_type']
    
    boundary_condition_left =  sim_perams['boundary_condition_left'] #incident_iso
    boundary_condition_right = sim_perams['boundary_condition_right'] #'reflecting'
    left_in_mag = sim_perams['left_in_mag']
    right_in_mag = sim_perams['right_in_mag']
    left_in_angle = sim_perams['left_in_angle']
    right_in_angle = sim_perams['right_in_angle']
    L = sim_perams['L']
    N_mesh = sim_perams['N_mesh']
    
    #snag some GL angles
    [angles_gq, weights_gq] = np.polynomial.legendre.leggauss(order_gauss_quad)
    
    #initilize some numpy soultion and problem space arrays
    angular_flux = np.zeros([order_gauss_quad, int(N_mesh*2)], data_type)
    angular_flux_next = np.zeros([order_gauss_quad, int(N_mesh*2)], data_type)
    angular_flux_last = np.zeros([order_gauss_quad, int(N_mesh*2)], data_type)
    
    scalar_flux  = np.zeros(int(N_mesh*2), data_type)
    scalar_flux_next = np.zeros(int(N_mesh*2), data_type)
    scalar_flux_last =  np.zeros(int(N_mesh*2), data_type)

    source_converged = False
    source_counter = 0
    no_convergence = False
    spec_rad = 0
    
    
    
    while source_converged == False:
        
        #print('Next Itteration: {0}'.format(source_counter),end='\r')
        
        #detemine bounds bounds for next itteration
        BCl = src.BoundaryCondition(boundary_condition_left,   0, N_mesh, angular_flux=angular_flux, incident_flux_mag=left_in_mag,  angle=left_in_angle,  angles=angles_gq)
        BCr = src.BoundaryCondition(boundary_condition_right, -1, N_mesh, angular_flux=angular_flux, incident_flux_mag=right_in_mag, angle=right_in_angle, angles=angles_gq)
        
        #find RHS of transport (see problem assigment)
        Q = src.RHSTransport(scalar_flux, xsec_scatter_mesh, source_mesh, N_mesh, dx_mesh)
        
        #simple corner balance to find angular flux for next itteration
        angular_flux = src.SCBRun(Q, xsec_mesh, dx_mesh, angles_gq, BCl, BCr, N_mesh)
        
        #calculate current
        current = src.Current(angular_flux, weights_gq, angles_gq)
        
        #calculate scalar flux for next itteration
        scalar_flux_next = src.ScalarFlux(angular_flux, weights_gq)
        
        if source_counter > 2:
            #check for convergence
            source_converged = src.HasItConverged(scalar_flux_next, scalar_flux)
            spec_rad = np.linalg.norm(scalar_flux_next - scalar_flux, ord=2) / np.linalg.norm((scalar_flux - scalar_flux_last), ord=2)
            
        #print()
        #print('Spec Rad {0}'.format(spec_rad))
        #print()
        
        
        #check for convergence
        source_converged = src.HasItConverged(scalar_flux_next, scalar_flux)
        
        #if stuck, display error then cut n run
        if source_counter > 10000:
            print('Error source not converged after 10000 itterations')
            print()
            source_converged = True
            no_convergence = True
        
        #reset for next itteration
        scalar_flux_last = scalar_flux
        scalar_flux = scalar_flux_next
        source_counter += 1
        
        
    #print()
    #print(source_counter)#spec_rad, no_convergence
    #print()scalar_flux, current
    
    if ret_angflux == True:
        scalar_flux = angular_flux_next
    
    return(scalar_flux, current, spec_rad, source_converged)





def OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, Q_tilde_modifier=None, time_dependent_mode=False):
    ''' Return converged scalar flux and current
    
    Implementation of simple corner balance source itterations
    '''
    
    #if running from command line print title contents
    #with open(os.path.join(sys.path[0], "title_print.txt"), "r", encoding="utf-8") as file:
    #    for line in file:
    #        print(line.strip())
    
    #Inputs: sets inputs from sim_perams tuple
    order_gauss_quad = sim_perams['N_angles']
    N_angles = order_gauss_quad
    half = int(N_angles/2)
    
    data_type = sim_perams['data_type']
    
    boundary_condition_left =  sim_perams['boundary_condition_left'] #incident_iso
    boundary_condition_right = sim_perams['boundary_condition_right'] #'reflecting'
    left_in_mag = sim_perams['left_in_mag']
    right_in_mag = sim_perams['right_in_mag']
    left_in_angle = sim_perams['left_in_angle']
    right_in_angle = sim_perams['right_in_angle']
    L = sim_perams['L']
    N_mesh = sim_perams['N_mesh']
    max_itter = sim_perams['max loops']
    
    #snag some GL angles
    [angles_gq, weights_gq] = np.polynomial.legendre.leggauss(order_gauss_quad)
    N = 2*N_mesh 
    
    #initilize some numpy soultion and problem space arrays
    angular_flux = np.zeros([order_gauss_quad, int(N)], data_type)
    angular_flux_next = np.zeros([order_gauss_quad, int(N)], data_type)
    angular_flux_last = np.zeros([order_gauss_quad, int(N)], data_type)
    
    scalar_flux  = np.zeros(int(N), data_type)
    scalar_flux_next  = np.ones(int(N), data_type)
    scalar_flux_last = np.ones(int(N), data_type)
    
    angular_flux_next = np.zeros([order_gauss_quad, int(N)], data_type)
    angular_flux_last = np.zeros([order_gauss_quad, int(N)], data_type)
    
    no_convergence = False
    source_converged = False
    source_counter = 0
    
    spec_rad_vec = np.empty([0])
    
    if time_dependent_mode==False:
        source_mesh = SourceMeshTransform(source_mesh, N_angles)
    
    assert (source_mesh.shape[0] == N_angles)
    assert (source_mesh.shape[1] == N)
    
    while source_converged == False:
        
        #print('Next Itteration: {0}'.format(source_counter),end='\r')
        
        #detemine bounds bounds for next itteration
        BCl = src.BoundaryCondition(boundary_condition_left,   0, N_mesh, angular_flux=angular_flux, incident_flux_mag=left_in_mag, angle=left_in_angle, angles=angles_gq)
        BCr = src.BoundaryCondition(boundary_condition_right, -1, N_mesh, angular_flux=angular_flux, incident_flux_mag=right_in_mag, angle=right_in_angle, angles=angles_gq)
        
        #simple corner balance to find angular flux for next itteration
        angular_flux_next = src.OCIRun(angular_flux, source_mesh, xsec_mesh, xsec_scatter_mesh, dx_mesh, angles_gq, weights_gq, BCl, BCr)
        
        #calculate current
        current = src.Current(angular_flux_next, weights_gq, angles_gq)
        
        #calculate scalar flux for next itteration
        scalar_flux_next = src.ScalarFlux(angular_flux_next, weights_gq)
        
        if source_counter > 2:
            #check for convergence
            source_converged = src.HasItConverged(angular_flux_next, angular_flux)
            spec_rad = np.linalg.norm(angular_flux_next - angular_flux, ord=2) / np.linalg.norm((angular_flux - angular_flux_last), ord=2)
        
        #if stuck, display error then cut n run
        if source_counter > max_itter:
            #print()
            #print('>>>WARNING: source not converged after 1000 itterations<<<')
            #print()
            no_convergence = True
            source_converged = True
        
        #reset for next itteration
        angular_flux_last = angular_flux
        angular_flux = angular_flux_next
        
        scalar_flux_last = scalar_flux
        scalar_flux = scalar_flux_next
        source_counter += 1
    
    #print(source_counter)
    #spec_rad, no_convergence
    
    if time_dependent_mode: #for time dependence
        scalar_flux = angular_flux_next
    
    return(scalar_flux, current, spec_rad, source_converged) #scalar_flux, current, 
    

def SourceMeshTransform(source_mesh_o, N_angles):
    '''Corrects the size of the source for transport when not time depenedent mode.
    Source must be [N_angles x 2*N_mesh] N_angles for nonisotropic, 2X for SCB
    '''
    
    N_mesh = source_mesh_o.size
    N_ans = int(2*N_mesh)
    source_mesh_n = np.zeros([N_angles, N_ans], np.float64)
    
    for i in range(N_mesh):
        for j in range(N_angles):
            source_mesh_n[j,2*i] = source_mesh_o[i]
            source_mesh_n[j,2*i+1] = source_mesh_o[i]
    
    return(source_mesh_n)

    
if __name__ == '__main__':
    x=0
    
    
    
    
