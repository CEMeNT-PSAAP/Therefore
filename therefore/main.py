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

def SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh):
    ''' Return converged scalar flux and current
    
    Implementation of simple corner balance source itterations
    '''
    
    #if running from command line print title contents
    with open(os.path.join(sys.path[0], "title_print.txt"), "r", encoding="utf-8") as file:
        for line in file:
            print(line.strip())
    
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
    
    scalar_flux  = np.zeros(int(N_mesh*2), data_type)
    scalar_flux_next  = np.zeros(int(N_mesh*2), data_type)

    source_converged = False
    source_counter = 0
    
    while source_converged == False:
        
        print('Next Itteration: {0}'.format(source_counter),end='\r')
        
        #detemine bounds bounds for next itteration
        BCl = src.BoundaryCondition(boundary_condition_left,   0, N_mesh, angular_flux=angular_flux, incident_flux_mag=left_in_mag, angle=left_in_angle, angles=angles_gq)
        BCr = src.BoundaryCondition(boundary_condition_right, -1, N_mesh, angular_flux=angular_flux, incident_flux_mag=right_in_mag, angle=right_in_angle, angles=angles_gq)
        
        #find RHS of transport (see problem assigment)
        Q = src.RHSTransport(scalar_flux, xsec_scatter_mesh, source_mesh, N_mesh, dx_mesh)
        
        #simple corner balance to find angular flux for next itteration
        angular_flux = src.SCBRun(angular_flux, Q, xsec_mesh, dx_mesh, angles_gq, BCl, BCr, N_mesh)
        
        #calculate current
        current = src.Current(angular_flux, weights_gq, angles_gq)
        
        #calculate scalar flux for next itteration
        scalar_flux_next = src.ScalarFlux(angular_flux, weights_gq)
        
        #check for convergence
        source_converged, error = src.HasItConverged(scalar_flux_next, scalar_flux)
        
        #if stuck, display error then cut n run
        if source_counter > 10000:
            print('Error source not converged after 1000 itterations')
            print()
            source_converged = True
        
        #reset for next itteration
        scalar_flux = scalar_flux_next
        source_counter += 1
    
    return(scalar_flux, current)
    
if __name__ == '__main__':
    x=0
