import numpy as np
import matplotlib.pyplot as plt

"Thomson's Rule for First-Time Telescope Makers: It is faster to make a four-inch mirror then a six-inch mirror than to make a six-inch mirror."

xsec = 10
scattering_ratio = .5
xsec_scattering = xsec*scattering_ratio

dx = 1
L = 3
N = int(L/dx)
N_mesh = int(2*N)
S = 0

dt = 0.1
N_time = 10
max_time = dt*(N_time-1)

#BCs incident iso
BCl = 10
BCr = 10

angular_flux      = np.zeros([2, N_mesh])
angular_flux_next = np.zeros([2, N_mesh])
angular_flux_midstep = np.zeros([2, N_mesh])
angular_flux_last = np.zeros([2, N_mesh])

angular_flux_final = np.zeros([2, int(N_mesh), N_time])

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

final_angular_flux_soultion = np.zeros([2, int(N_mesh)])

for k in range(N_time):

    angular_flux      = np.zeros([2, N_mesh])
    angular_flux_last = np.zeros([2, N_mesh])
    angular_flux_midstep = np.zeros([2, N_mesh])
    angular_flux_midstep_last = np.zeros([2, N_mesh])

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

            A = np.array([[-mu1/2 - w1*manaz + gamma, mu1/2,                    -w2*manaz,                 0],
                        [-mu1/2,                    -mu1/2 - w1*manaz + gamma,  0,                       -w2*manaz],
                        [-w1*manaz,                 0,                         mu2/2 + gamma - w2*manaz, mu2/2],
                        [0,                         -w1*manaz,                 -mu2/2,                   mu2/2 + gamma - w2*manaz]])

            if i == 0: #left bc
                b = np.array([[dx/4*S],
                            [dx/4*S - mu1 * angular_flux[0, i*2+2]],
                            [dx/4*S + mu2 * BCl],
                            [dx/4*S]])
            elif i == N-1: #right bc
                b = np.array([[dx/4*S],
                            [dx/4*S - mu1 * BCr],
                            [dx/4*S + mu2 * angular_flux[1, i*2-1]],
                            [dx/4*S]])
            else: #mid communication
                b = np.array([[dx/4*S],
                            [dx/4*S - mu1 * angular_flux[0, i*2+2]],
                            [dx/4*S + mu2 * angular_flux[1, i*2-1]],
                            [dx/4*S]])
            
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
plt.show()
#plt.savefig('Test Angular flux')

#
def A_neg():
    gamma = (dx*xsec_total_time)/2
    timer = dx/(v*dt)
    timer2 = dx/(2*v*dt)
    a = mu/2

    np.array([[-a + gamma, a,          timer2,            0],
              [-a,         -a + gamma, 0,                 timer2],
              [-timer2,    0,          timer - a + gamma, a],
              [0,          -timer,     -a,                timer -a + gamma]])

def A_pos():
    gamma = (dx*xsec_total_time)/2
    timer = dx/(v*dt)
    timer2 = dx/(2*v*dt)
    a = mu/2

    np.array([[a + gamma, a, timer2, 0],
              [-a, a + gamma, 0, timer2],
              [-timer, 0, timer + a + gamma, a],
              [0, -timer, -a, timer +a + gamma]])

def c_neg():
    timer2 = dx/(2*v*dt)

    np.array([[dx/4],
              [],
              [],
              []])

def c_pos():
    timer2 = dx/(2*v*dt)

    np.array([[dx/4*Ql + timer2*psi_halfLast_L + mu * psi_leftBound],
              [dx/4*Qr + timer2*psi_halfLast_R],
              [dx/4*Q_halfNext_L + mu*psi_halfNext_leftBound],
              [dx/4*Q_halfNext_R - mu*psi_halfNext_rightBound]])

def scatter_source():

    np.array([[],
              [],
              [],
              []])