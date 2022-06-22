"""
Created on Sun May 15 20:34:03 2022

@author: jacksonmorgan
"""

import numpy as np
import matplotlib.pyplot as plt
import therefore
from timeit import default_timer as timer

def flatLinePlot(x, y):
    for i in range(y.size):
        xx = x[i:i+2]
        yy = [y[i], y[i]]
        plt.plot(xx, yy, '-k')

data_type = np.float64

L = 10
xsec = 10


#dx = np.linspace(.01, 2, 75)
dx = .5

ratio = np.linspace(0, .99, 75)

#source = 0
source = np.linspace(0, 10, 75)

mfp = xsec*dx
#print(mfp)

xs = ratio.size
ys = source.size


no_converge_oci = np.zeros([xs, ys])
spec_rad_oci = np.zeros([xs, ys])
no_converge_si = np.zeros([xs, ys])
spec_rad_si = np.zeros([xs, ys])

total_runs = xs * ys

time_oci = 0
time_si = 0

for i in range(xs):
    for k in range(ys):
        
        print('Percent done: %2d' %(((i*ys+k)/total_runs)*100),end='\r')
        
        scattering_xsec = xsec*ratio[i]
        
        N_mesh = int(L/dx)

        in_flux = source[k]/((xsec-scattering_xsec)/2)

        dx_mesh = dx*np.ones(N_mesh, data_type)
        xsec_mesh = xsec*np.ones(N_mesh, data_type)
        xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
        source_mesh = source[k]*np.ones(N_mesh, data_type)

        sim_perams = {'data_type': data_type,
                      'N_angles': 2,
                      'L': L,
                      'N_mesh': N_mesh,
                      'boundary_condition_left': 'incident_iso',
                      'boundary_condition_right': 'incident_iso',
                      'left_in_mag': 10,
                      'right_in_mag': 10,
                      'left_in_angle': .3,
                      'right_in_angle': -.3,
                      'max loops': 10000}

        start = timer()
        [sf, cur, spec_rad_oci[i,k], no_converge_oci[i,k]] = therefore.OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
        time_oci += timer() - start
        start = timer()
        [sf, cur, spec_rad_si[i,k], no_converge_si[i,k]] = therefore.SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
        time_si += timer() - start
        #print(spec_rad[i,k])

if ((no_converge_si == True).any):
    print()
    print('>>>WARNING: Some runs of SI did not converge before itter kick.')
    print('            Recomend grow max_itter value and run again  <<<')
    print('')

if ((no_converge_oci == True).any):
    print()
    print('>>>WARNING: Some runs of OCI did not converge before itter kick.')
    print('            Recomend grow max_itter value and run again  <<<')
    print('')

print()
print('Total time for OCI computations: {0}'.format(time_oci))
print('Total time for SI computations:  {0}'.format(time_si))
print()

mfp = xsec*dx
#np.savez('spectral_radius_s4', mfp=mfp, ratio=ratio, spec_rad_oci=spec_rad_oci, spec_rad_si=spec_rad_si)

[Xx, Yy] = np.meshgrid(source,ratio)


fig = plt.figure(1)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Xx,Yy,spec_rad_oci)
plt.title('OCI SCB Spectral Radius Plot')
plt.xlabel('mfp [σ*Δx]')
plt.ylabel('Scattering Ratio [$σ_s$/σ]')
ax.set_zlabel('Spectrial Radius [ρ]')
plt.show()


fig = plt.figure(2)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Xx,Yy,spec_rad_si)
plt.title('SI SCB Spectral Radius Plot')
plt.xlabel('mfp [σ*Δx]')
plt.ylabel('Scattering Ratio [$σ_s$/σ]')
ax.set_zlabel('Spectrial Radius [ρ]')
plt.show()

