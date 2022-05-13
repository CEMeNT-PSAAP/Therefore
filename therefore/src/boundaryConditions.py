import numpy as np

def BoundaryCondition(BC, i, N_mesh, angular_flux=None, incident_flux_mag=None, angles=None):

    if BC == 'reflecting':
        N_angles: int = angular_flux.shape[0]
        half = int(N_angles/2)
        psi_required = np.zeros(N_angles)
        
        print(half)
        
        if i == 0: #left edge
            psi_required[half:] = angular_flux[:half, 0]
            
        else:
            psi_required[:half] = angular_flux[half:, -1]

        
    elif BC == 'vacuum':
        psi_required = np.zeros(angular_flux.shape[0])
    
    elif BC == 'incident_iso':
        psi_required = BC_isotropic(incident_flux_mag)
        
    elif BC == 'incident_ani':
        psi_required = BC_ani(incident_flux_mag, angle)
        
    else:
        print()
        print('>>>Error: No Boundary Condition Supplied<<<')
        print()
    
    return(psi_required)

def BC_isotropic(incident_flux_mag):
    BC = (incident_flux_mag/angles_gq.size)*np.ones(angles_gq.size)
    return(BC)

def BC_ani(incident_flux_mag, angle):
    angle_id = find_nearest(angles_gq, angle)
    BC = np.zeros(angles_gq.size)
    BC[angle_id] = incident_flux_mag
    return(BC)
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return (idx)

def RHS_transport(scalar_flux, scattering_xsec, source, N_mesh, dx):
    Q = np.zeros(N_mesh*2)
    for i in range(N_mesh):
        Q[2*i]   = scalar_flux[2*i]   * scattering_xsec[i]/2 + (source/2)*dx/2
        Q[2*i+1] = scalar_flux[2*i+1] * scattering_xsec[i]/2 + (source/2)*dx/2
    return(Q)
