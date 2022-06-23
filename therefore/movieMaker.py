import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def MovieMaker(scalar_flux, L, FLP=False):
    '''You like Qinton Ternintino? Damn, we didn't measure this in feet
    '''
    
    N_mesh = scalar_flux.shape[0]
    N_time = scalar_flux.shape[1]
    
    t_gif = 5 #how long the gif should lasts
    frames_per = int(N_time/t_gif)
    
    
    t = np.linspace(0, 10, N_time+1)
    x = np.linspace(0, L, int(N_mesh))
    
    fig,ax = plt.subplots() #plt.figure(figsize=(6,4))
    
    ax.grid()
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'Flux')
    ax.set_title(r'$\bar{\phi}_{k,j}$')
    
    line1, = ax.plot(x, scalar_flux[:,0],'-k',label="Therefore")
    #line2, = ax.plot([], [],'--r',label="AZURV1")
    text   = ax.text(0.02, 0.9, '', transform=ax.transAxes)
    ax.legend()  
          
    def animate(k):
        line1.set_ydata(scalar_flux[:,k])
        #line2.set_data(my_xx, my_sf[:,k])
        text.set_text(r'$t \in [%.1f,%.1f]$ s'%(t[k],t[k+1]))
        return line1,#, text
    
    simulation = animation.FuncAnimation(fig, animate, frames=N_time)
    plt.show()
    
    writervideo = animation.PillowWriter(fps=frames_per)
    simulation.save('slab_reactor.gif') #saveit!


#ax.fill_between(x_mid,phi[k,:]-phi_sd[k,:],phi[k,:]+phi_sd[k,:],alpha=0.2,color='b')
#text.set_text(r'$t \in [%.1f,%.1f]$ s'%(t[k],t[k+1]))
