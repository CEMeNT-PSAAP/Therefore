import numpy as np
import matplotlib.pyplot as plt

xsec = 10
scattering_ratio = .5
xsec_scattering = xsec*scattering_ratio

dx = 1
L = 3
N = int(L/dx)
N_mesh = 2*N
S = 0

#BCs incident iso
BCl = 10
BCr = 10

#print('Fuck')

angular_flux      = np.zeros([2, int(N_mesh)])
angular_flux_next = np.zeros([2, int(N_mesh)])
angular_flux_last = np.zeros([2, int(N_mesh)])

mu1 = -0.57735
mu2 = 0.57735

w1 = 1
w2 = 1

tol = 1e-6
error = 1
max_itter = 1000
itter = 0

manaz = dx*xsec_scattering/4
gamma = xsec*dx/2

while error > tol or max_itter < itter:

    print()
    print("========================================")
    print("next cycle")
    print("========================================")
    print()

    # TODO: OCI
    for i in range(N):
        A = np.zeros([4,4])
        b = np.zeros([4,1])

        A = np.array([[-mu1/2 - w1*manaz, -mu1/2 + gamma,    -w2*manaz,                0],
                      [-mu1/2 + gamma,    mu1/2 - w1*manaz,  0,                        -w2*manaz],
                      [-w1*manaz,         0,                 mu2/2 + gamma - w2*manaz, mu2/2],
                      [0,                 -w1*manaz,         -mu2/2,                   mu2/2 + gamma - w2*manaz]])

        if i == 0: #left bc
            b = np.array([[dx/2*S - mu1 * angular_flux[0, i*2+2]],
                          [dx/2*S],
                          [dx/2*S - mu2 * BCl],
                          [dx/2*S]])
        elif i == N-1: #right bc
            b = np.array([[dx/2*S - mu1 * BCr],
                          [dx/2*S],
                          [dx/2*S - mu2 * angular_flux[1, i*2-1]],
                          [dx/2*S]])
        else: #mid communication
            b = np.array([[dx/2*S - mu1 * angular_flux[0, i*2+2]],
                          [dx/2*S],
                          [dx/2*S - mu2 * angular_flux[1, i*2-1]],
                          [dx/2*S]])
        
        print("Large cell %d".format(i))
        print(b)
        print()
        print(A)
        print()

        angular_flux_next[:,2*i:2*i+2] = np.linalg.solve(A,b).reshape(-1,2)
        
        print(angular_flux_next[:,2*i:2*i+2])
        print()

        itter += 1 

    # TODO: Error
    error = np.linalg.norm(angular_flux_next - angular_flux, ord=2)

    angular_flux_last = angular_flux
    angular_flux = angular_flux_next
    
f=1
X = np.linspace(0, L, int(N_mesh))
plt.figure(f)
plt.plot(X, angular_flux[0,:],  '-*k',  label='OCI 1')
plt.plot(X, angular_flux[1,:],  '--*k', label='OCI 2')
#plt.plot(X, scalar_flux2[0,:], '-r',  label='SI 1')
#plt.plot(X, scalar_flux2[1,:], '--r', label='SI 2')
plt.title('Test Flux')
plt.xlabel('Distance')
plt.ylabel('Angular Flux')
plt.savefig('Test Angular flux')

#
