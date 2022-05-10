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
    return(np.sum(weights * angular_flux))

def Current(angular_flux, angles)
    return(angular_flux * angles)

def HasItConverged(scalar_flux_next, scalar_flux):
   np.allclose(scalar_flux_next, scalar_flux, rtol=tol)


SimpleCorn


#simple corner balance
def SimpleCornerBalanceRight(Q_l, Q_r, psi_mh, xsec, dx, mu):
    
    laguz  = (xsec * dx) / 4
    mannaz = mu/2 - laguz
    othala = laguz + mu/2
    
    denominator = mu/2 + (mannaz*mu)/(2*othala) + laguz*(mannaz/othala + 1)
    
    psi_l = (dx/2*Q_l - ((mu*dx)/(4*othala) + (xsec*dx**2)/(6*othala))*Q_r + mu*psi_mh) / denominator
    
    psi_r = psi_l*mannaz/othala + dx/(2*othala)*Q_r
    
    return(psi_l, psi_r)
    
def SimpleCornerBalanceLeft():
    Q_l = 
    Q_r = 
    
    laguz  = (xsec * dx) / 4
    mannaz = mu/2 - laguz
    othala = laguz + mu/2
    
    denominator = mu/2 + (mannaz*mu)/(2*othala) + laguz*(mannaz/othala + 1)
    
    psi_l = (dx/2*Q_l - ((mu*dx)/(4*othala) + (xsec*dx**2)/(6*othala))*Q_r + mu*psi_mh) / denominator
    
    psi_r = psi_ml*mannaz/othala + dx/(2*othala)*Q_r
    
    return(psi_ml, psi_mr)

def ExitingFlux():
    if mu > 0:
        if i == 0:
            psi = L_bound
        elif i > 0:
            psi = 

#print title contents
with open("title_print.txt", "r", encoding="utf-8") as file:
    for line in file:
        print(line.strip())


#Inputs

order_gauss_quad = 4
scattering_xsec = 
total_xsec = 
dx = 
incident_scalar_flux_left
incident_scalar_flux_right = 
boundary_condition_left = 0
boundary_condition_right = 0


# TODO: Mesh building



[weights_gq, angles_gq] = np.polynomial.legendre.leggauss(order_gauss_quad)

scalar_flux = np.zeros(N_cells)

source_convergence = False

# TODO: Source itterations
while source_convergence == False:
    
    RHS_transport = (scalar_flux * scattering_xsec)/2 + source/2
    
    # TODO: Discrete Ordinants: Sweep Left to Right (+)
    for mu in range(N_angles):
        for i in range(N_cells):
            angular_flux_next[mu, i] = SimpleCornerBalanceRight(Q_l, Q_r, psi_mh, xsec, dx, mu)
            
    
    # TODO: Discrete Ordinants: Sweep Right to Left (-)  
    for mu in range(N_angles)
        for i in range(N_cells):
            angular_flux_next[mu, i] = SimpleCornerBalance()
    
    scalar_flux_next = ScalarFlux(angular_flux_next, weights_gq)
    
    # TODO: Check for convergence
    source_convergence = HasItConverged(scalar_flux_next, scalar_flux)
    
# TODO: Negativie flux fixups

# TODO: Plot scalar flux and current
