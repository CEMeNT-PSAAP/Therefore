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
        #print('i in scalar flux comp {0}'.format(i))
        scalar_flux += weights[i] * angular_flux[i,:]
        
    return(scalar_flux)

def Current(angular_flux, angles, weights):
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
            for i in range(N_mesh-1, -1, -1):
                #check bound
                if i == N_mesh-1:
                    psi_mh = BCr[angle]
                else:
                    psi_mh = psi_l[i+1]
                
                # inputs switch around when going to the left
                # in kernel l to r and r to l
                
                #entering  exiting               entering  exiting
                [psi_r[i], psi_l[i]] = SCBKernel_Linalg(Q[2*i+1], Q[2*i], psi_mh, xsec[i], dx[i], mu[angle])
                '''
                print('--------------------------------------------')
                print('Angle: {0}   Cell: {1}'.format(mu[angle], i))
                print()
                print('Q_l:    {0}'.format(Q[2*i+1]))
                print('Q_r:    {0}'.format(Q[2*i]))
                print('psi_mh: {0}'.format(psi_mh))
                print('xsec:   {0}'.format(xsec[i]))
                print('dx:     {0}'.format(dx[i]))
                print('mu:     {0}'.format(mu[angle]))
                print('psi_l:  {0}'.format(psi_l[i]))
                print('psi_r:  {0}'.format(psi_r[i]))
                print()
                print()
                '''
        else:
            for i in range(N_mesh):
                #print(i)
                #check bound
                if i == 0:
                    psi_mh = BCl[angle]
                else:
                    psi_mh = psi_r[i-1]
                
                [psi_l[i], psi_r[i]] = SCBKernel_Linalg(Q[2*i], Q[2*i+1], psi_mh, xsec[i], dx[i], mu[angle])
                '''
                print('--------------------------------------------')
                print('Angle: {0}   Cell: {1}'.format(mu[angle], i))
                print()
                print('Q_l:    {0}'.format(Q[2*i+1]))
                print('Q_r:    {0}'.format(Q[2*i]))
                print('psi_mh: {0}'.format(psi_mh))
                print('xsec:   {0}'.format(xsec[i]))
                print('dx:     {0}'.format(dx[i]))
                print('mu:     {0}'.format(mu[angle]))
                print('psi_l:  {0}'.format(psi_l[i]))
                print('psi_r:  {0}'.format(psi_r[i]))
                print()
                print()
                '''
        #print('Angle! {0}'.format(mu[angle]))
        #print()
        ##print(psi_l)
        #print()
        #print(psi_r)
        #print()
        
        for i in range(N_mesh):
            angular_flux[angle, 2*i]   = psi_l[i]
            angular_flux[angle, 2*i+1] = psi_r[i]
        
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
    
    b = np.array([[dx/2 * Q_exiting],[dx/2*Q_entering - mu*psi_mh]])
    
    [psi_exiting, psi_entering] = np.linalg.solve(A,b)
    
    return(psi_entering, psi_exiting)

def BoundaryCondition(BC, i, N_mesh, angular_flux=None, incident_flux_mag=None, angle=None, angles=None):

    if BC == 'reflecting':
        N_angles: int = angular_flux.shape[0]
        half = int(N_angles/2)
        psi_required = np.zeros(N_angles)
        
        if i == 0:
            psi_required[:half] = angular_flux[half:, -1]
        else:
            psi_required[half:] = angular_flux[:half, 0]

        
    elif BC == 'vacuum':
        psi_required = np.zeros(angular_flux.shape[0])
    
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
    dx = 0.01
    L = 1
    N_mesh = int(L/dx)
    incident_flux_mag = [0,0]
    incidnet_flux_angle = [0,0]
    boundary_condition_right = 'reflecting' #'reflecting'
    boundary_condition_left =  'reflecting'
    source = 1
    
    
    data_type = np.float64

    # TODO: Mesh building
    dx_mesh = dx*np.ones(N_mesh, data_type)
    xsex_mesh = xsec*np.ones(N_mesh, data_type)
    scattering_xsec = scattering_xsec*np.ones(N_mesh, data_type)
    np.ones(N_mesh, data_type)
    
    [angles_gq, weights_gq] = np.polynomial.legendre.leggauss(order_gauss_quad)
    
    # angles_gq = np.cos(angles_gq)
    
    #print(weights_gq)
    #print(angles_gq)

    angular_flux = np.zeros([order_gauss_quad, int(N_mesh*2)], data_type)
    scalar_flux  = np.zeros(int(N_mesh*2), data_type)
    scalar_flux_next  = np.zeros(int(N_mesh*2), data_type)


    # TODO: Source itterations
    source_converged = False
    source_counter = 0
    BCr = np.zeros(2)
    BCl = np.zeros(2)
    
    while source_converged == False:
    
        print('Next Itteration! {0}'.format(source_counter))
        #print('Angular flux {0}'.format(angular_flux))
        #print()
        #print()
        #BCr = angular_flux[:,-1]
        #BCl = angular_flux[:,0]
        
        #BCr = BoundaryCondition(boundary_condition_left, -1, N_mesh, angular_flux, angles=angles_gq)
        #print('BCr {0}'.format(BCr))
        
        #BCl = BoundaryCondition(boundary_condition_right, 0, N_mesh, angular_flux, angles=angles_gq)
        #print('BCl {0}'.format(BCl))
        
        BCr = angular_flux[:, -1]
        #BCr[1] = angular_flux[0, -1]
        
        BCl = angular_flux[:, 0]
        #BCl[1] = angular_flux[0, 0]
        
        Q = RHS_transport(scalar_flux, scattering_xsec, source, N_mesh, dx)
        
        # TODO: simple corner balance
        angular_flux = SCBRun(angular_flux, Q, xsex_mesh, dx_mesh, angles_gq, BCl, BCr, N_mesh)
        
        # TODO: calculate current
        current = Current(angular_flux, weights_gq, angles_gq)
        
        # TODO: calculate scalar flux for next itteration
        scalar_flux_next = ScalarFlux(angular_flux, weights_gq)
        
        
        # TODO: Check for convergence
        source_converged, error = HasItConverged(scalar_flux_next, scalar_flux)
        
        if source_counter > 1000:
            print('Error source not converged after 1000 itterations')
            print()
            source_converged = True
        
        #print()
        #print(error)
        #print()
        #print(scalar_flux_next)
        """
        
        print()
        print()
        print('================================================')
        print()
        print()
        print(scalar_flux_next)
        """
        
        scalar_flux = scalar_flux_next
        source_counter += 1
    
    print()
    print('Source counter: {0}'.format(source_counter))
    print()
    
    # TODO: Negativie flux fixups
    # not required for balance methods

    #print(source_counter)
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
