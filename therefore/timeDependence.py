import numpy as np
from .itterationSchemes import SourceItteration, OCI
import therefore.src as src

def TimeLoop(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source, theta=1):
    velocity = sim_perams['velocity']
    dt = sim_perams['dt']
    N_time = sim_perams['N_time']
    N_angles = sim_perams['N_angles']
    N_mesh = sim_perams['N_mesh']
    data_type = sim_perams['data_type']
    
    N_ans = int(2*N_mesh)
    angular_flux_total = np.zeros([N_angles, N_ans, N_time], data_type)
    angular_flux_last = np.zeros([N_angles, N_ans], data_type)
    current_total = np.zeros([N_ans, N_time], data_type)
    scalar_flux = np.zeros([N_ans, N_time], data_type)
    spec_rad = np.zeros(N_time)
    
    [angles, weights] = np.polynomial.legendre.leggauss(N_angles)
    
    source_mesh = np.ones([N_angles, N_ans], data_type)
    for i in range(N_mesh):
        for j in range(N_angles):
            source_mesh[j,i] = source[i]
            source_mesh[j,i+1] = source[i]
    
    
    for t in range(N_time):
        xsec_mesh_t = xsec_mesh + (1/(velocity* theta* dt))
        source_mesh_tilde = source_mesh + angular_flux_last/(velocity* theta* dt)
        
        [angular_flux_total[:,:,t], current_total[:,t], spec_rad[t], source_converged] = OCI(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh_tilde, angular_flux_last, True)
        
        if source_converged == False:
            print()
            print('>>>WARNING<<<')
            print('   Method of itteration did not converge!')
            print('')
            
        #angular_flux_total[:,:,t] = TimeDiscretization(angular_flux_last, angular_flux_half)
        
        print('Time loop! {0}'.format(t))
        
        angular_flux_last = angular_flux_total[:,:,t]
        scalar_flux[:,t] = src.ScalarFlux(angular_flux_last, weights)
    
    return(scalar_flux, current_total, spec_rad)
    
    
    
    
#def TimeDiscretization(angular_flux_last, angular_flux_half, theta):
#    '''Using diamond discretization in time
#    '''
    
#    angular_flux_next = (angular_flux_half - theta*angular_flux_last) / (1-theta)
    
#    return(angular_flux_next)
