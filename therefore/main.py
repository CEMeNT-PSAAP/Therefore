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

#utility functions
def ScalarFlux(angular_flux, weights):
    scalar_flux = np.zeros(angular_flux.shape[1])
    for i in range(angular_flux.shape[0]):
        print('i in scalar flux comp {0}'.format(i))
        scalar_flux += weights[i] * angular_flux[i,:]
        
    return(scalar_flux)

def Current(angular_flux, angles, weights):a
    current = np.zeros(angular_flux.shape[1])
    for i in range(angular_flux.shape[0]):
        current += weights[i] * angles[i] * angular_flux[i,:]
    return(current)

def HasItConverged(scalar_flux_next, scalar_flux, tol=1e-8):
   error = max(abs((scalar_flux_next - scalar_flux) / scalar_flux_next))
   return(error < tol, error)

#Simple Corner balence sweep
def SCBRun(angular_flux, Q, xsec, dx, mu, BCl, BCr, N_mesh):
    
    for angle in range(mu.size):
        psi_l = np.zeros(N_mesh, np.float64)
        psi_r = np.zeros(N_mesh, np.float64)
    
        if mu[angle] < 0:
            for i in range(N_mesh-1, 0, -1):
                #check bound
                if i == N_mesh-1:
                    psi_mh = BCr[angle]
                else:
                    psi_mh = psi_l[i+1]
                
                # inputs switch around when going to the left
                # in kernel l to r and r to l
                [psi_r[i], psi_l[i]] = SCBKernel(Q[2*i+1], Q[2*i], psi_mh, xsec[i], dx[i], mu[angle])
        
        else:
            for i in range(N_mesh):
                print(i)
                #check bound
                if i == 0:
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
    
    mannaz = mu/2 + xsec*dx/2
    
    denominator = mannaz + mu/(4*mannaz)
    
    psi_l = (dx/2*Q_l - (mu*dx)/(4*mannaz)*Q_r + mu*psi_mh) / denominator
    
    psi_r = (Q_r*dx/2 + psi_l/2) / mannaz
    
    return(psi_l, psi_r)


def BoundaryCondition(BC, i, N_mesh, angular_flux=None, incident_flux_mag=None, angle=None, angles=None):

    if BC == 'reflecting':
        N_angles: int = angular_flux.size[0]
        half: int = N_angles/2
        psi_required = np.zeros(N_angles)
        
        for j in range(N_angles.size[0])
            if angles[j] > 0:
                psi_required[0, j] = angular_flux[i, j + half]
            else:
                psi_required[-1, j] = angular_flux[-1, j - half]
        
    elif BC == 'vacuum':
        psi_required = np.zeros(angular_flux.size[0])
    
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

def RHS_transport(scalar_flux, scattering_xsec, source, N_mesh, dx):
    Q = np.zeros(N_mesh*2)
    for i in range(N_mesh):
        Q[2*i]   = scalar_flux[2*i]   * scattering_xsec[i]/2 + (source/2)*dx/2
        Q[2*i+1] = scalar_flux[2*i+1] * scattering_xsec[i]/2 + (source/2)*dx/2
    return(Q)


def main():
    #print title contents
    with open("title_print.txt", "r", encoding="utf-8") as file:
        for line in file:
            print(line.strip())


    #Inputs

    order_gauss_quad = 2
    scattering_xsec = 0
    xsec = 3
    dx = 0.25
    L = 1
    N_mesh = int(L/dx)
    incident_flux_mag = [0,0]
    incidnet_flux_angle = [0,0]
    boundary_condition_right = 'reflecting'
    boundary_condition_left =  'reflecting'
    source = 1
    
    
    data_type = np.float64

    # TODO: Mesh building
    dx_mesh = dx*np.ones(N_mesh, data_type)
    xsex_mesh = xsec*np.ones(N_mesh, data_type)
    scattering_xsec = scattering_xsec*np.ones(N_mesh, data_type)
    np.ones(N_mesh, data_type)
    
    [angles_gq, weights_gq] = np.polynomial.legendre.leggauss(order_gauss_quad)
    
    print(weights_gq)
    print(angles_gq)

    angular_flux = np.zeros([order_gauss_quad, int(N_mesh*2)], data_type)
    scalar_flux  = np.zeros(int(N_mesh*2), data_type)
    scalar_flux_next  = np.zeros(int(N_mesh*2), data_type)


    # TODO: Source itterations
    source_converged = False
    source_counter = 0
    while source_converged == False:
        print('Angular flux {0}'.format(angular_flux))
        print()
        print()
        #BCr = angular_flux[:,-1]
        #BCl = angular_flux[:,0]
        
        BCr = BoundaryCondition(boundary_condition_left, -1, N_mesh, angular_flux, angles=angles_gq)
        print('BCr {0}'.format(BCr))
        
        BCl = BoundaryCondition(boundary_condition_right, 0, N_mesh, angular_flux, angles=angles_gq)
        print('BCl {0}'.format(BCl))
        
        Q = RHS_transport(scalar_flux, scattering_xsec, source, N_mesh, dx)
        print('Q {0}'.format(Q))
        print()
        print()
        
        # TODO: simple corner balance
        angular_flux = SCBRun(angular_flux, Q, xsex_mesh, dx_mesh, angles_gq, BCl, BCr, N_mesh)
        
        # TODO: calculate current
        current = Current(angular_flux, weights_gq, angles_gq)
        
        # TODO: calculate scalar flux for next itteration
        scalar_flux_next = ScalarFlux(angular_flux, weights_gq)
        
        
        # TODO: Check for convergence
        source_converged, error = HasItConverged(scalar_flux_next, scalar_flux)
        
        print()
        print(error)
        print()
        print(scalar_flux_next)
        print()
        print()
        print(scalar_flux)
        
        
        scalar_flux = scalar_flux_next
        source_counter += 1
        
    # TODO: Negativie flux fixups
    # not required for balance methods

    print(source_counter)
    print()

    # TODO: Plot scalar flux and current
    print('plotting')
    X = np.arange(N_mesh*2)
    plt.figure(1)
    plt.plot(X, scalar_flux)
    plt.show()
    print('plotted')
    
if __name__ == '__main__':
    main()
