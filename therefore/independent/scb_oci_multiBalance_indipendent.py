import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=np.inf)

"""Thomson's Rule for First-Time Telescope Makers: It is faster to make a 
four-inch mirror then a six-inch mirror than to make a six-inch mirror."""

def A_neg(dx, v, dt, mu, xsec_total):
    gamma = (dx*xsec_total)/2
    timer = dx/(v*dt)
    timer2 = dx/(2*v*dt)
    a = mu/2

    A_n = np.array([[-a + gamma, a,          timer2,            0],
                    [-a,         -a + gamma, 0,                 timer2],
                    [-timer,     0,          timer - a + gamma, a],
                    [0,          -timer,     -a,                timer -a + gamma]])
    
    return(A_n)



def A_pos(dx, v, dt, mu, xsec_total):
    gamma = (dx*xsec_total)/2
    timer = dx/(v*dt)
    timer2 = dx/(2*v*dt)
    a = mu/2

    A_p = np.array([[a + gamma, a,         timer2,            0],
                    [-a,        a + gamma, 0,                 timer2],
                    [-timer,    0,         timer + a + gamma, a],
                    [0,         -timer,    -a,                timer +a + gamma]])

    return(A_p)



def c_neg(dx, v, dt, mu, Ql, Qr, Q_halfNext_L, Q_halfNext_R, psi_halfLast_L, psi_halfLast_R, psi_rightBound, psi_halfNext_rightBound):
    timer2 = dx/(2*v*dt)

    c_n = np.array([[dx/4*Ql + timer2*psi_halfLast_L],
                    [dx/4*Qr + timer2*psi_halfLast_R - mu* psi_rightBound],
                    [dx/4*Q_halfNext_L],
                    [dx/4*Q_halfNext_R - mu*psi_halfNext_rightBound]])

    return(c_n)



def c_pos(dx, v, dt, mu, Ql, Qr, Q_halfNext_L, Q_halfNext_R, psi_halfLast_L, psi_halfLast_R, psi_leftBound, psi_halfNext_leftBound):
    timer2 = dx/(2*v*dt)

    c_p = np.array([[dx/4*Ql + timer2*psi_halfLast_L + mu * psi_leftBound],
                    [dx/4*Qr + timer2*psi_halfLast_R],
                    [dx/4*Q_halfNext_L + mu*psi_halfNext_leftBound],
                    [dx/4*Q_halfNext_R]])
    
    return(c_p)



def scatter_source(dx, xsec_scattering, N, w):
    S = np.zeros([2*N,2*N])
    beta = dx*xsec_scattering/4
    for i in range(N):
        for j in range(N):
            S[2*i,2*j] = beta * w[i]
            S[2*i+1,2*j+1] = beta * w[i]

    return(S)

xsec = 10
scattering_ratio = 0.0
xsec_scattering = xsec*scattering_ratio

printer = False
printer_TS = True

dx = .1
L = 1
N = int(L/dx)
N_mesh = int(2*N)
Q = 0

dt = 0.25
N_time = 5
max_time = dt*(N_time-1)
v = 1

#BCs incident iso
BCl = 10
BCr = 0

angular_flux      = np.zeros([2, N_mesh])
angular_flux_next = np.zeros([2, N_mesh])
angular_flux_midstep = np.zeros([2, N_mesh])
angular_flux_last = np.zeros([2, N_mesh])

angular_flux_final = np.zeros([2, int(N_mesh), N_time])

mu1 = -0.57735
mu2 =  0.57735

w1 = 1
w2 = 1
w = np.array([w1, w2])

N_angle = 2

tol = 1e-9
error = 1
max_itter = 1000

manaz = dx*xsec_scattering/4
gamma = xsec*dx/2

final_angular_flux_solution = np.zeros([N_time, N_angle, N_mesh])
final_angular_flux_midstep_solution = np.zeros([N_time, N_angle, N_mesh])



# the zeroth stored solution is the initial condition
for k in range(1, N_time, 1):

    if (printer_TS):
        print()
        print("========================================")
        print("next time step!")
        print("========================================")
        print()

    # iterating on these till convergence
    angular_flux      = np.zeros([2, N_mesh]) 
    angular_flux_last = np.zeros([2, N_mesh])   # last refers to last iteration
    angular_flux_midstep = np.zeros([2, N_mesh])
    angular_flux_midstep_last = np.zeros([2, N_mesh])   # last refers to last iteration

    #initial guesses?
    itter = 0
    error_eos = 1
    while error_eos > tol and max_itter > itter:
        print(itter)
        if itter == max_itter:
            print('Crap: {0}'.format(k))

        if (printer):
            print()
            print("========================================")
            print("next cycle: {0}".format(itter))
            print("========================================")
            print()

        # TODO: OCI
        for i in range(N):
            i_l = int(2*i)
            i_r = int(2*i+1)

            A = np.zeros([8,8])
            c = np.zeros([8,1])
            
            A_small = A_neg(dx, v, dt, mu1, xsec)
            S_small = scatter_source(dx, xsec_scattering, N_angle, w)
            #print('>>> S_small <<<')
            #print(S_small)
            #print()

            assert ((A_small.size == S_small.size))

            A[:4, :4] = A_small - S_small

            A_small = A_pos(dx, v, dt, mu2, xsec)
            S_small = scatter_source(dx, xsec_scattering, N_angle, w)

            assert ((A_small.size == S_small.size))

            A[4:, 4:] = A_small - S_small

            psi_halfLast_L = final_angular_flux_midstep_solution[k-1, :, i_l] # known
            psi_halfLast_R = final_angular_flux_midstep_solution[k-1, :, i_r] # known

            # boundary conditions
            if i == 0:  #left
                psi_rightBound = angular_flux_last[0, i_r+1] # iterating on (unknown)
                psi_leftBound =  BCl # known

                psi_halfNext_rightBound = angular_flux_midstep_last[0, i_r+1] # iterating on (unknown)
                psi_halfNext_leftBound  = BCl # known

            elif i == N-1: #right
                psi_rightBound = BCr # known
                psi_leftBound =  angular_flux_last[1, i_l-1] # iterating on (unknown)

                psi_halfNext_rightBound = BCr # known
                psi_halfNext_leftBound  = angular_flux_midstep_last[1, i_l-1] # iterating on (unknown)

            else: #middles
                psi_rightBound = angular_flux_last[0, i_r+1] # iterating on (unknown)
                psi_leftBound =  angular_flux_last[1, i_l-1] # iterating on (unknown)

                psi_halfNext_rightBound = angular_flux_midstep_last[0, i_r+1] # iterating on (unknown)
                psi_halfNext_leftBound  = angular_flux_midstep_last[1, i_l-1] # iterating on (unknown)


            c[:4] = c_neg(dx, v, dt, mu1, Q, Q, Q, Q, psi_halfLast_L[0], psi_halfLast_R[0], psi_rightBound, psi_halfNext_rightBound)
            c[4:] = c_pos(dx, v, dt, mu2, Q, Q, Q, Q, psi_halfLast_L[1], psi_halfLast_R[1], psi_leftBound, psi_halfNext_leftBound)

            if (printer):
                print("Large cell {0}".format(i))
                print('>>> psi_right bound (BC) <<<')
                print(psi_rightBound)
                print()
                print('>>> c vector <<<')
                print(c)
                print()
                print('>>> full A mat <<<')
                print(A)
                print()

            angular_flux_raw = np.linalg.solve(A,c)

            # resorting into proper locations in solution vectors
            for p in range(N_angle):
                angular_flux[p,i_l]         = angular_flux_raw[4*p]
                angular_flux[p,i_r]         = angular_flux_raw[4*p+1]
                
                angular_flux_midstep[p,i_l] = angular_flux_raw[4*p+2]
                angular_flux_midstep[p,i_r] = angular_flux_raw[4*p+3]

            if (printer):
                print('>>> raw solution <<<')
                print(angular_flux_raw)
                print()
                print('>>> angular flux eos reorganized <<<')
                print(angular_flux)
                print()
                print('>>> angular flux mid step reorganized <<<')
                print(angular_flux_midstep)
                print()


            itter += 1 

        # TODO: Error
        error_eos = np.linalg.norm(angular_flux_midstep - angular_flux_midstep_last, ord=2)
        error_mos = np.linalg.norm(angular_flux - angular_flux_last, ord=2)

        final_angular_flux_solution[k, :, :] = angular_flux
        final_angular_flux_midstep_solution[k, :, :] = angular_flux_midstep
            
        angular_flux_last = angular_flux 
        angular_flux_midstep_last = angular_flux_midstep
        
        

f=1
X = np.linspace(0, L, int(N_mesh))
plt.figure(f)
plt.plot(X, final_angular_flux_solution[1, 1,:],  '--*g', label='0')
#plt.plot(X, final_angular_flux_midstep_solution[1, 1,:],  '-*g',  label='0 + 1/2')
plt.plot(X, final_angular_flux_solution[2, 1,:],  '--*k', label='1')
#plt.plot(X, final_angular_flux_midstep_solution[2, 1,:],  '-*k',  label='1 + 1/2')
plt.plot(X, final_angular_flux_solution[3, 1,:],  '--*r', label='2')
#plt.plot(X, final_angular_flux_midstep_solution[3, 1,:],  '-*r',  label='2 + 1/2')
plt.plot(X, final_angular_flux_solution[4, 1,:],  '--*b', label='9, 2')
#plt.plot(X, final_angular_flux_midstep_solution[-1, 1,:],  '-*b',  label='3 + 1/2')
#plt.plot(X, scalar_flux2[0,:], '-r',  label='SI 1')
#plt.plot(X, scalar_flux2[1,:], '--r', label='SI 2')
plt.title('Test Ang Flux: Positive ordinant')
plt.xlabel('Distance')
plt.ylabel('Angular Flux')
plt.legend()
plt.show()
#plt.savefig('Test Angular flux')