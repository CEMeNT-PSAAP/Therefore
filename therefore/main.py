"""
Therefore
Method of implementation for slab wall problem
Author: Jackson P. Morgan
breif: An implementaiton of simplified corner balance for assignemtn 2 in NSE 653
prof: Todd Palmer
date: May 9th 2022
"""


import numpy as np
import matplotlib as pyplot

#utility functions
def ScalarFlux(angular_flux, weights):
    return(np.sum(weights * angular_flux, axis=1))

def Current(angular_flux, angles)
    return(np.sum(weights * angular_flux * angles, axis=1))

def HasItConverged(scalar_flux_next, scalar_flux):
   np.allclose(scalar_flux_next, scalar_flux, rtol=tol)

#Simple Corner balence sweep
def SCBRun(angular_flux, Q, xsec, dx, mu, BCl, BCr):
    
    for angle in range(mu.size):
    psi_l = np.zeros(N_mesh, data_type)
    psi_r = np.zeros(N_mesh, data_type)
    
        if angle < 0:
            for i in range(-N_mesh):
            
                #check bound
                if i == N_mesh-1:
                    psi_mh = BCr[angle]
                else:
                    psi_mh = psi_l[i+1]
                
                [psi_r[i], psi_l[i]] = SCBKernel(Q[2*i+1], Q[2*i], psi_mh, xsec[i], dx[i], mu[angle])
        
        else:
            for i in range(N_mesh):
                
                #check bound
                if i == N_mesh-1:
                    psi_mh = BCl[angle]
                else:
                    psi_mh = psi_r[i-1]
                
                [psi_l[i], psi_r[i]] = SCBKernel(Q[2*i], Q[2*i+1], psi_mh, xsec[i], dx[i], mu[angle])
                
        for i in range(N_mesh):
            angular_flux[angle, 2*i]   = psi_l[i]
            angular_flux[angle, 2*i+1] = psi_r[i]
        
    return(angular_flux)
    

#simple corner balance for a single cell
def SCBKernel(Q_l, Q_r, psi_mh, xsec, dx, mu):
    
    laguz  = (xsec * dx) / 4
    mannaz = mu/2 - laguz
    othala = laguz + mu/2
    
    denominator = mu/2 + (mannaz*mu)/(2*othala) + laguz*(mannaz/othala + 1)
    
    psi_l = (dx/2*Q_l - ((mu*dx)/(4*othala) + (xsec*dx**2)/(6*othala))*Q_r + mu*psi_mh) / denominator
    
    psi_r = psi_l*mannaz/othala + dx/(2*othala)*Q_r
    
    return(psi_l, psi_r)


def BoundaryCondition(BC, i):

    if BC == 'reflecting':
        if i == 1:
            i = N_mesh-1
        psi_required = angular_flux[:,i]
        
    elif BC == 'incident_iso':
        psi_required = BC_isotropic(incident_flux_mag[i])
        
    elif BC == 'incident_ani':
        psi_required = BC_ani(incident_flux_mag[i], angle[i])
    else:
        print()
        print('>>>Error: No Boundary Condition Supplied<<<')
        print()
    
    return(psi_required)

def BC_isotropic(incident_flux_mag):
    BC = (incident_flux_mag/angles_gq.size)*np.ones(angles_gq.size)
    return(BC)

def BC_ani(incident_flux_mag, angle):
    angle_id = find_nearest(angles_gq, angle)
    BC = np.zeros(angles_gq.size)
    BC[angle_id] = incident_flux_mag
    return(BC)
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return (idx)

def RHS_transport():
    Q = np.zeros(N_mesh*2)
    for i in range(N_mesh):
        Q[2*i]   = scalar_flux[2*i]   * scattering_xsec[i]/2 + source/2
        Q[2*i+1] = scalar_flux[2*i+1] * scattering_xsec[i]/2 + source/2
    return(Q)


def main():
    #print title contents
    with open("title_print.txt", "r", encoding="utf-8") as file:
        for line in file:
            print(line.strip())


    #Inputs

    order_gauss_quad = 4
    scattering_xsec = 1
    xsec = 1
    dx = 0.1
    L = 1
    N_mesh = int(L/dx)
    incident_flux_mag = [0,0]
    incidnet_flux_angle = [0,0]
    boundary_condition_right = 'reflecting'
    boundary_condition_left =  'reflecting'
    source = 1
    
    
    data_type = np.float64

    # TODO: Mesh building
    dx = dx*np.ones(N_mesh, data_type)
    xsex = = xsec*np.ones(N_mesh, data_type)
    scattering_xsec = scattering_xsec*np.ones(N_mesh, data_type)
    np.ones(N_mesh, data_type)
    
    [weights_gq, angles_gq] = np.polynomial.legendre.leggauss(order_gauss_quad)


    angular_flux = np.zeros([order_gauss_quad, int(N_mesh*2)], data_type)
    scalar_flux  = np.zeros(int(N_mesh*2), data_type)
    scalar_flux_next  = np.zeros(int(N_mesh*2), data_type)


    # TODO: Source itterations
    source_converged = False
    while source_converged == False:
        
        BCr = BoundaryCondition(boundary_condition_left, 0)
        BCl = BoundaryCondition(boundary_condition_right, 1)
        
        Q = RHS_transport()
        
        # TODO: simple corner balance
        angular_flux = SCBRun(scalar_flux, angular_flux, Q, xsec, dx, angles_gq, BCl, BCr)
        
        # TODO: calculate current
        current = Current(angular_flux_next, weights_gq, angles_gq)
        
        # TODO: calculate scalar flux for next itteration
        scalar_flux_next = ScalarFlux(angular_flux_next, weights_gq)
        
        # TODO: Check for convergence
        source_converged = HasItConverged(scalar_flux_next, scalar_flux)
        
        scalar_flux = scalar_flux_next
        
    # TODO: Negativie flux fixups
    # not required for balance methods


    # TODO: Plot scalar flux and current
    
if __name__ == '__main__':
    main()
