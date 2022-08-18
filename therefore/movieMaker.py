import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from .azurv1 import azurv1_spav

def MovieMaker(scalar_flux, scalar_flux2, L, FLP=False):
    '''You like Qinton Ternintino? Damn, we didn't measure this in feet
    '''
    
    N_mesh = scalar_flux.shape[0]
    N_time = scalar_flux.shape[1]
    
    t_gif = 5 #how long the gif should lasts
    frames_per = int(N_time/t_gif)
    
    
    t = np.linspace(0, 10, N_time+1)
    x = np.linspace(0, L, int(N_mesh))
    x_eval = np.linspace(-L/2, L/2, int(N_mesh+1))
    
    fig,ax = plt.subplots() #plt.figure(figsize=(6,4))
    
    ax.grid()
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'Flux')
    ax.set_title(r'$\bar{\phi}_{k,j}$')
    
    line1, = ax.plot(x, scalar_flux[:,0], '-k',label="Therefore SI")
    line2, = ax.plot(x, scalar_flux2[:,0],'-r',label="Therefore OCI")
    line3, = ax.plot(x, azurv1_spav(x_eval, .9, 0.01), '--b',label="AZURV1")
    text   = ax.text(0.02, 0.9, '', transform=ax.transAxes)
    ax.legend()  
    plt.ylim([0,2])
    plt.xlim([4.9, 5.1])
    
    def animate(k):
        line1.set_ydata(scalar_flux[:,k])
        line2.set_ydata(scalar_flux2[:,k])
        line3.set_ydata(azurv1_spav(x_eval, .9, .1*k+.01))
        text.set_text(r'$t \in [%.1f,%.1f]$ s'%(t[k],t[k+1]))
        print('Figure production percent done: {0}'.format(int(k/N_time)*100), end = "\r")
        return line2, line2, line3, #, text
    print()
    print()
    
    simulation = animation.FuncAnimation(fig, animate, frames=N_time)
    plt.show()
    
    writervideo = animation.PillowWriter(fps=frames_per)
    simulation.save('slab_reactor.gif') #saveit!


#ax.fill_between(x_mid,phi[k,:]-phi_sd[k,:],phi[k,:]+phi_sd[k,:],alpha=0.2,color='b')
#text.set_text(r'$t \in [%.1f,%.1f]$ s'%(t[k],t[k+1]))
