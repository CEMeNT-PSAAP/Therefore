import numpy as np
import numba as nb
import scipy.integrate as inter



#thanks Aaron!
def azurv1(x, c, t):
    '''Produce the anaylitical soultion for a vector
    returns a vector
    '''
    N_mesh = x.size
    scalar_flux = np.zeros(N_mesh)
    
    for i in range(N_mesh):
        scalar_flux[i] = phi(x[i], c, t)
    
    return(scalar_flux)
    


#time average
def avurv1_TimeSpaceAvg(xGrid, tGrid, c):
    #assuming constant spacining
    dx = xGrid[1] - xGrid[0]
    dt = tGrid[1] - tGrid[0]

    Nx = xGrid.size-1
    Nt = tGrid.size-1

    scalar_flux_avg = np.zeros([Nx, Nt])
    scalar_flux_grid = np.zeros([Nx, Nt+1])

    print(tGrid[0]-dt)
    print(tGrid[-1]+dt)
    print(Nt+1)
    tGrid_avg = np.linspace(tGrid[0]-dt, tGrid[-1]+dt, Nt+1)

    print(tGrid_avg)
    print(xGrid)
    print()
    print(c)
    print()
    half = int(Nx/2)

    #true average in space
    for i in range(Nt):
        print(tGrid_avg[i])
        scalar_flux_grid[:,i] = azurv1_spav(xGrid, tGrid_avg[i], c)
        print(i)
        print(scalar_flux_grid[:,i])
        print(scalar_flux_grid[half,i])
        print(i)

    #simple average in time
    #for i in range(Nt):
    #    for j in range(Nx):
    #        scalar_flux_avg[j,i] = (scalar_flux_grid[j,i] + scalar_flux_grid[j,i+1]) / 2

    assert(scalar_flux_avg.shape[0] == Nx)
    assert(scalar_flux_avg.shape[1] == Nt)

    return(scalar_flux_avg)


def azurv1_spav(x, c, t):
    '''Spatially averaged! x vector is location of cell bounds (so N_mesh+1)
    '''
    N_mesh = x.size-1
    dx = x[1] - x[0] #assuming spatially independednt
    scalar_flux = np.zeros(N_mesh)
    
    for i in range(N_mesh):
        scalar_flux[i] = inter.quad(phi, x[i], x[i+1], args=(c, t), limit=1000)[0] / dx
    
    return(scalar_flux)



# Evaluate flux at position x and time t fir a scatter ratio c
def phi(x, c, t):
    if t == 0.0 or abs(x) >= t:
        return 0.0
        
    eta = x / t
    
    if c > 0:
        integral = inter.quad(integrand,0,np.pi,args=(eta,t,c),limit=1000)[0]
    else:
        integral = 0
        
    return (np.exp(-t) / (2 * t)) *(1 + (c * t / (4 * np.pi)) * (1 - eta ** 2) * integral)



# Helper function for evaluating flux
def integrand(u,eta,t,c):
    i = complex(0,1)
    
    q   = (1+eta)/(1-eta)
    
    xi = (np.log(q)+i*u)/(eta+i*np.tan(u/2))
    return (1.0 / (np.cos(u / 2))) ** 2 * (xi ** 2 * np.exp(((c * t) / 2) * (1 - eta ** 2) * xi)).real
    
    

if __name__ == "__main__":

    x = np.arange(0, 1, .01)
    t = np.arange(0, 1, .02)
    c = .9

    scalar_flux_avg = avurv1_TimeSpaceAvg(x, t, c)

    print(scalar_flux_avg.shape)