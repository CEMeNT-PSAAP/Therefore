"""
Created on Sun May 15 20:34:03 2022
@author: jacksonmorgan
"""

import numpy as np
import matplotlib.pyplot as plt
import therefore


data_type = np.float64

L = 1
dx = 0.5
xsec = 1
ratio = 0
scattering_xsec = xsec*ratio
source_mat = 0
source_a = 2
N_mesh = int(L/dx)

dx_mesh = dx*np.ones(N_mesh, data_type)
xsec_mesh = xsec*np.ones(N_mesh, data_type)
xsec_scatter_mesh = scattering_xsec*np.ones(N_mesh, data_type)
source_mesh = source_mat*np.ones(N_mesh, data_type)

#psi_in = source_mat / (xsec*(1-ratio)/2)

sim_perams = {'data_type': data_type,
              'N_angles': 2,
              'L': L,
              'N_mesh': N_mesh,
              'boundary_condition_left':  'incident_iso',
              'boundary_condition_right': 'incident_iso',
              'left_in_mag': 10,
              'right_in_mag': 10,
              'left_in_angle': .3,
              'right_in_angle': 1,
              'max loops': 10000,
              'tolerance': 1e-9}

[scalar_flux, current, spec_rad, conver, loops] = therefore.OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
[scalar_flux2, current2, spec_rad2, conver2, loops] = therefore.SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)