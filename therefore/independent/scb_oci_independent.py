import numpy as np
import matplotlib.pyplot as plt




xsec =
scattering_ratio = 
xsec_scattering = 

dx = 1
L = 3
N = int(L/dx)
N_mesh = 2*N
S = 0

#BCs incident iso
BCl = 10
BCr = 10


angular_flux      = np.zeros([2, int(N_mesh)], data_type)
angular_flux_next = np.zeros([2, int(N_mesh)], data_type)
angular_flux_last = np.zeros([2, int(N_mesh)], data_type)

mu1 = -0.57735
mu2 = 0.57735


tol = 1e-6
error = 1
max_itter = 1000
itter = 0

manaz = dx*xsec_scattering/4
gamma = xsec*dx/2

while error > tol and max_itter < itter:

    # TODO: OCI
    for i in range(N):
        A = np.zeros([4,4])
        b = np.zeros([4,1])

        A = np.array([[-mu1/2 - w1*manaz, -mu1/2 + gamma, -w2*manaz, 0],
                      [-mu1/2 + gamma, mu1/2 - w1*manaz, 0, -w2*manaz],
                      [-w1*manaz, 0, mu2/2 + gamma - w2*manaz, mu2/2],
                      [0, -w1*manaz, -mu2/2, mu2/2 + gamma - w2*manaz]])

        if i = 0: #left bc
            b = np.array([[dx/2*S - mu1 * 0],
                          [dx/2*S],
                          [dx/2*S - m2 * BCl],
                          [dx/2*S]])
        elif i = N-1: #right bc
            b = np.array([[dx/2*S - m1 * BCr],
                          [dx/2*S],
                          [dx/2*S - m2 * 0],
                          [dx/2*S]])
        else: #mid communication
            b = np.array([[dx/2*S - m1 * angular_flux[]],
                          [dx/2*S],
                          [dx/2*S - m2 * angular_flux[]],
                          [dx/2*S]])
        

    # TODO: Error


    angular_flux_last = angular_flux
    angular_flux = angular_flux_next
        

#
