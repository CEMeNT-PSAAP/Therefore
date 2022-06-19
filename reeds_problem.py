"""
Created on Sun May 15 20:34:03 2022

@author: jacksonmorgan
"""

import numpy as np
import matplotlib.pyplot as plt
import therefore


def flatLinePlot(x, y, pl):
    for i in range(y.size):
        xx = x[i:i+2]
        yy = [y[i], y[i]]
        plt.plot(xx, yy, pl)

data_type = np.float64

#mesh builder for Reed's Problem
region_id = np.array([1,2,3,4,5], int)
region_widths = np.array([2,2,1,1,2], int)
region_bounds = np.array([2,4,5,6,8], float)
sigma_s = np.array([.9, .9, 0, 0, 0], data_type)
sigma_t = np.array([1, 1, 0, 5, 50], data_type)
Source = np.array([0, 1, 0, 0, 50], data_type)
dx = .1/sigma_t
dx[2] = .25  #fix a nan
N_region = np.array(region_widths/dx, int)

N_mesh: int = sum(N_region)

xsec_mesh = np.empty(N_mesh, data_type)
xsec_scatter_mesh = np.empty(N_mesh, data_type)
dx_mesh = np.empty(N_mesh, data_type)
source_mesh = np.empty(N_mesh, data_type)
region_id_mesh = np.empty(N_mesh, data_type)
region_id_mesh_2 = np.empty(N_mesh*2, data_type)

#build the mesh
for i in range(region_widths.size):
    LB = sum(N_region[:i])
    RB = sum(N_region[:i+1])
    xsec_mesh[LB:RB] = sigma_t[i]
    xsec_scatter_mesh[LB:RB] = sigma_s[i]
    dx_mesh[LB:RB] = dx[i]
    source_mesh[LB:RB] = Source[i]
    region_id_mesh[LB:RB] = region_id[i]

for i in range(N_mesh):
    region_id_mesh_2[2*i] = region_id_mesh[i]
    region_id_mesh_2[2*i+1] = region_id_mesh[i]
    
#set sim peramweters
sim_perams = {'data_type': data_type,
              'N_angles': 8,
              'L': 0,
              'N_mesh': N_mesh,
              'boundary_condition_left': 'vacuum',
              'boundary_condition_right': 'reflecting',
              'left_in_mag': 0,
              'right_in_mag': 0,
              'left_in_angle': 0,
              'right_in_angle': 0,
              'max loops': 10000}

#launch source itterations
[scalar_flux, current, spec_rad, source_converged] = therefore.SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
[scalar_flux2, current2, spec_rad2, source_converged] = therefore.OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)

print()
print('Overall spectral radius of SI: {0}'.format(spec_rad))
print('Overall spectral radius of OCI: {0}'.format(spec_rad2))
print()

#post process and plot
Y_reg = [0, max(scalar_flux)*2]
X_reg = [0,0]

f = 1
plt.figure(f)

x_plot = np.zeros(N_mesh*2+1)


for i in range(N_mesh):
    x_plot[2*i] = sum(dx_mesh[:i])
    x_plot[2*i+1] = sum(dx_mesh[:i+1])

x_plot[-1] = 8.00001

#scalar_flux2 /= 2

flatLinePlot(x_plot, scalar_flux, '-k')
flatLinePlot(x_plot, scalar_flux2, '-r')
#plt.title('Problem 2: Reeds Problem')
plt.xlabel('Distance [cm]')
plt.ylabel('Scalar Flux')
plt.ylim([0,max(scalar_flux)*1.25])

# plot region demarkers
for i in range(5):
    X_reg[0] = region_bounds[i]
    X_reg[1] = region_bounds[i]
    plt.plot(X_reg, Y_reg, c='lightgrey')
plt.ylim([0,8])
plt.show()


f = 2
plt.figure(f)


flatLinePlot(x_plot, current, '-k')
flatLinePlot(x_plot, current2, '-r')
#plt.title('Problem 2: Reeds Problem')
plt.xlabel('Distance [cm]')
plt.ylabel('Current')
plt.ylim([0,max(scalar_flux)*1.25])

# plot region demarkers
for i in range(5):
    X_reg[0] = region_bounds[i]
    X_reg[1] = region_bounds[i]
    plt.plot(X_reg, Y_reg, c='lightgrey')
plt.ylim([0,3])
plt.show()

