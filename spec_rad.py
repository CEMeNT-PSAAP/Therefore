"""
Created on Sun May 15 20:34:03 2022

@author: jacksonmorgan
"""

import numpy as np
import matplotlib.pyplot as plt
import therefore


def flatLinePlot(x, y):
    for i in range(y.size):
        xx = x[i:i+2]
        yy = [y[i], y[i]]
        plt.plot(xx, yy, '-k')

data_type = np.float64

L = 10
xsec = 10
source = 0

dx = np.linspace(.01, 2, 75)
ratio = np.linspace(0, .99, 75)

mfp = xsec*dx
#print(mfp)

no_converge = np.zeros([ratio.size, dx.size])
spec_rad = np.zeros([ratio.size, dx.size])

total_runs = ratio.size * dx.size

for i in range(ratio.size):
    for k in range(dx.size):
        
        print('Percent done: %2d' %(((i*dx.size+k)/total_runs)*100),end='\r')
        
        scattering_xsec = xsec*ratio[i]
        
        N_mesh = int(L/dx[k])

        in_flux = source/((xsec-scattering_xsec)/2)

        dx_mesh = dx[k]*np.ones(N_mesh, data_type)
        xsec_mesh = xsec*np.ones(N_mesh, data_type)
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

        #launch source itterations #SourceItteration
        [spec_rad[i,k], no_converge[i,k]] = therefore.OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
        #print(spec_rad[i,k])

if ((no_converge == True).any):
    print()
    print('>>>WARNING: Some runs did not converge before itter kick.')
    print('            Recomend grow max_itter value and run again  <<<')
    print('')

mfp = xsec*dx
np.savez('spectral_radius_s4', mfp=mfp, ratio=ratio, spec_rad=spec_rad)

[Xx, Yy] = np.meshgrid(mfp,ratio)

fig = plt.figure(1)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Xx,Yy,spec_rad)
plt.title('SI SCB Spectral Radius Plot')
plt.xlabel('mfp [σ*Δx]')
plt.ylabel('Scattering Ratio [$σ_s$/σ]')
ax.set_zlabel('Spectrial Radius [ρ]')
plt.show()

