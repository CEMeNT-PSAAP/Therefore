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

v = 2

#dx = np.linspace(.01, 2, 75)
dx = .5

ratio = np.linspace(0, .99, 45)

source = 0
#mfp = np.linspace(0, 20, 75)

dt = np.linspace(0.0001,1,45)

#mfp = xsec*dx
#print(mfp)

xs = ratio.size
ys = dt.size


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
        
        xsec_hat = xsec + 1/(v*dt[k])
        scattering_xsec = xsec*ratio[i]
        
        N_mesh = int(L/dx)

        #print(scattering_xsec)
        #print(dx)
        #print(N_mesh)
        #in_flux = source[k]/((xsec-scattering_xsec)/2)

        dx_mesh = dx*np.ones(N_mesh, data_type)
        xsec_mesh = xsec_hat*np.ones(N_mesh, data_type)
        xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
        source_mesh = source*np.ones(N_mesh, data_type)

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
        [sf, cur, spec_rad_oci[i,k], no_converge_oci[i,k], loops] = therefore.OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
        time_oci += timer() - start
        start = timer()
        [sf, cur, spec_rad_si[i,k], no_converge_si[i,k], loops] = therefore.SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
        time_si += timer() - start
        #print(spec_rad[i,k])

if ((no_converge_si == False).any):
    np.set_printoptions(linewidth=np.inf)
    print()
    print('>>>WARNING: Some runs of SI did not converge before itter kick.')
    print('            Recomend grow max_itter value and run again  <<<')
    print('')
    print(no_converge_si)
    print()

if ((no_converge_oci == False).any):
    np.set_printoptions(linewidth=np.inf)
    print()
    print('>>>WARNING: Some runs of OCI did not converge before itter kick.')
    print('            Recomend grow max_itter value and run again  <<<')
    print('')
    print(no_converge_oci)
    print()

print()
print('Total time for OCI computations: {0}'.format(time_oci))
print('Total time for SI computations:  {0}'.format(time_si))
print()


mfp = xsec*dx
#mfp = xsec*dx
np.savez('spectral_radius_s2', mfp=mfp, ratio=ratio, spec_rad_oci=spec_rad_oci, spec_rad_si=spec_rad_si)

[Xx, Yy] = np.meshgrid(dt,ratio)



fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(Xx,Yy,spec_rad_oci)
plt.title('OCI SCB Spectral Radius Plot')
plt.xlabel('dt [σ*Δx]')
plt.ylabel('Scattering Ratio [$σ_s$/σ]')
ax.set_zlabel('Spectrial Radius [ρ]')
plt.savefig('specrad_oci',dpi=600)


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(Xx,Yy,spec_rad_si)
plt.title('SI SCB Spectral Radius Plot')
plt.xlabel('dt')
plt.ylabel('Scattering Ratio [$σ_s$/σ]')
ax.set_zlabel('Spectrial Radius [ρ]')
plt.savefig('specrad_si',dpi=600)

epsilon = 1e-9
N_si = np.zeros([xs,ys])
N_oci = np.zeros([xs,ys])

for i in range(xs):
    for k in range(ys):
        N_si[i,k] = np.log(epsilon)/np.log(spec_rad_si[i,k])
        N_oci[i,k] = np.log(epsilon)/np.log(spec_rad_oci[i,k])

N_ratio = N_si/N_oci

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(Xx,Yy,N_ratio)
plt.title('Ratio')
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'Scattering Ratio [$\sigma_s/\sigma$]')
ax.set_zlabel(r'$N_{SI}/N_{OCI}$')

for ii in range(0,360,1):
        ax.view_init(elev=10., azim=ii)
        plt.savefig("figs/movie%d.png" % ii)

plt.savefig('ratio_value',dpi=600)