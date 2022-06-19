import numpy as np
from .itterationSchemes import SourceItteration, OCI


def TimeLoop(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh, theta=1):
    velocity = sim_perams['velocity']
    dt = sim_perams['dt']
    N_time = sim_perams['max time']
    N_angles = sim_perams['N_angles']
    N_mesh = sim_perams['N_mesh']
    
    angular_flux_total = np.null([N_angles, N_mesh, N_time])
    current_total = np.null([N_mesh, N_time])
    spec_rad = np.null(N_time)
    
    for t in range(N_time):
        xsec_mesh_t = xsec_mesh + (1/(velocity* theta* dt))
        source_mesh_t = source_mesh + angular_flux_last/(velocity* theta* dt)
        [angular_flux_total[:,:,t], current_total[:,t], spec_rad[t], source_converged] = therefore.SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
        
        if source_converged == False:
            print()
            print('>>>WARNING<<<')
            print('   Method of itteration did not converge!')
            print('')
        
    return(angular_flux_total, current_total, spec_rad)
    
    
    
    
def TimeDiscretization(angular_flux_last, angular_flux_half, theta):
    '''Using diamond discretization in time
    '''
    
    angular_flux_next = (angular_flux_half - theta*angular_flux_last) / (1-theta)
    
    return(angular_flux_next)
