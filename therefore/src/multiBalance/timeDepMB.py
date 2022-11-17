import numpy as np

from timeit import default_timer as timer
from .scb_oci_mb import OCIMBTimeStep
from .scb_si_mb import OCIMBSITimeStep
import therefore.src.utilities as utl


def multiBalance(inital_angular_flux, sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source, backend='OCI_MB'):
    velocity = sim_perams['velocity']
    dt = sim_perams['dt']
    N_time = sim_perams['N_time']
    N_angles = sim_perams['N_angles']
    N_mesh = sim_perams['N_mesh']
    data_type = sim_perams['data_type']
    printer = sim_perams['print']
    
    N_ans = int(2*N_mesh)
    current_total = np.zeros([N_ans, N_time], data_type)
    scalar_flux = np.zeros([N_ans, N_time+1], data_type)
    spec_rad = np.zeros(N_time)

    [angles, weights] = np.polynomial.legendre.leggauss(N_angles)

    #sotring the intial scalar_flux
    scalar_flux[:,0] = utl.ScalarFlux(inital_angular_flux, weights)
    
    source_mesh = np.ones([N_angles, N_ans], data_type)
    for i in range(N_mesh):
        for j in range(N_angles):
            source_mesh[j,2*i] = source[i]
            source_mesh[j,2*i+1] = source[i]
    


    angular_flux_mid_last = inital_angular_flux
    angular_flux = np.zeros([N_angles, int(2*N_mesh), N_time])

    angular_flux_mid = np.zeros([N_angles, int(2*N_mesh), N_time])

    source_converged: bool = True
    
    for t in range(N_time):
        
        start = timer()

        if (backend == 'OCI_MB'):
            [angular_flux[:,:,t], angular_flux_mid[:,:,t], current_total[:,t], spec_rad[t], loops, 
            source_converged] = OCIMBTimeStep(sim_perams, angular_flux_mid_last, source_mesh, xsec_mesh, xsec_scatter_mesh, dx_mesh, angles, weights)
        elif (backend == 'SI_MB'):
            [angular_flux[:,:,t], angular_flux_mid[:,:,t], current_total[:,t], spec_rad[t], 
            loops, source_converged] = OCIMBSITimeStep(sim_perams, angular_flux_mid_last, source_mesh, xsec_mesh, xsec_scatter_mesh, dx_mesh, angles, weights)
        else:
            print('>>>ERROR: NO Backend provided')
            print('     select between OCI and SI!')
            print()
            
        end = timer()
        
        
        if source_converged == False:
            print()
            print('>>>WARNING<<<')
            print('   Method of iteration did not converge!')
            print('')
            
        #angular_flux_total[:,:,t] = TimeDiscretization(angular_flux_last, angular_flux_half)
        
        #psi_check = source_mesh_tilde[0,5] / (xsec_mesh_t[5]*(1)/2)
        
        if (printer):
            print('Time step: {0}'.format(t))
            print('     -ρ:          {0}    '.format(spec_rad[t]))
            print('     -wall time:  {0} [s]'.format(end-start))
            print('     -loops:      {0}'.format(loops))
            print()
        #print('     -Ψ check:    {0}    '.format(psi_check))
        
        
        angular_flux_mid_last = angular_flux_mid[:,:,t]
        scalar_flux[:,t+1] = utl.ScalarFlux(angular_flux_mid_last, weights)
    
    return(scalar_flux, current_total, spec_rad)
    
    
#def TimeDiscretization(angular_flux_last, angular_flux_half, theta):
#    '''Using diamond discretization in time
#    '''
    
#    angular_flux_next = (angular_flux_half - theta*angular_flux_last) / (1-theta)
    
#    return(angular_flux_next)
