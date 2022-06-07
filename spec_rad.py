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

dx = np.linspace(.1, 1.1, 50)
ratio = np.linspace(0, .75, 50)

#print(dx*xsec)

converged = np.zeros([ratio.size, dx.size])
spec_rad = np.zeros([ratio.size, dx.size])

total_runs = ratio.size * dx.size

for i in range(ratio.size):
    for k in range(dx.size):
        
        print('Percent done: {0}'.format(((i*ratio.size+k)/total_runs)*100),end='\r')
        
        scattering_xsec = xsec*ratio[i]
        
        N_mesh = int(L/dx[k])

        in_flux = source/((xsec-scattering_xsec)/2)

        #print()
        #print('Soultion should be of magnitude: {0}'.format(in_flux))
        #print()

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
                      'right_in_angle': -.3}

        #launch source itterations #SourceItteration
        [spec_rad[i,k], converged[i,k]] = therefore.SourceItteration(sim_perams, dx_mesh, xsec_mesh, xsec_scatter_mesh, source_mesh)
print()
print()

print(spec_rad)
print()
print()
print()
print()
print(converged)





mfp = xsec*dx
[Xx, Yy] = np.meshgrid(mfp,ratio)

fig = plt.figure(1)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Xx,Yy,spec_rad)
plt.title('OCI SCB Spectral Radius Plot')
plt.xlabel('mfp [σ*Δx]')
plt.ylabel('Scattering Ratio [$σ_s$/σ]')
ax.set_zlabel('Spectrial Radius [ρ]')

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#ax.xaxis.set_major_locator(LinearLocator(10))
#ax.xaxis.set_major_formatterplt.clabel(surf,fontsize=10)#,fontcolor='k')(FormatStrFormatter('%2.0f'))
plt.show()

np.savez('spectral_radius_s4', mfp=mfp, ratio=ratio, spec_rad=spec_rad)

'''
data = [[ 66386, 174296,  75131, 577908,  32015],
        [ 58230, 381139,  78045,  99308, 160454],
        [ 89135,  80552, 152558, 497981, 603535],
        [ 78415,  81858, 150656, 193263,  69638],
        [139361, 331509, 343164, 781380,  52269]]
columns = ('Freeze', 'Wind', 'Flood', 'Quake', 'Hail')
rows = ['%d year' % x for x in (100, 50, 20, 10, 5)]
values = np.arange(0, 2500, 500)
value_increment = 1000
# Get some pastel shades for the colors
colors = plt.cm.BuPu(np.linspace(0, 0.5, len(mfp)))
n_rows = len(data)
index = np.arange(len(columns)) + 0.3
bar_width = 0.4
# Initialize the vertical-offset for the stacked bar chart.
y_offset = np.zeros(len(columns))
# Plot bars and create text labels for the table
cell_text = []
for row in range(mfp.size):
    plt.bar(index, spec_rad[row], bar_width, bottom=y_offset, color=colors[row])
    y_offset = y_offset + spec_rad[row]
    cell_text.append(['%1.1f' % (x / 1000.0) for x in y_offset])
# Reverse colors and text labels to display the last value at the top.
colors = colors[::-1]
cell_text.reverse()
# Add a table at the bottom of the axes
the_table = plt.table(cellText=cell_text,
                      rowLabels=mfp,
                      rowColours=colors,
                      colLabels=columns,
                      loc='bottom')
# Adjust layout to make room for the table:
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.ylabel("Loss in ${0}'s".format(value_increment))
plt.yticks(values * value_increment, ['%d' % val for val in values])
plt.xticks([])
plt.title('Loss by Disaster')
# plt.show()
'''





'''
f=1
X = np.linspace(0, L, int(N_mesh))
plt.figure(f)
plt.plot(X, scalar_flux)
plt.title('Infinte Med')
plt.xlabel('Distance')
plt.ylabel('Scalar Flux')
plt.show()

f=1
X = np.linspace(0, L, int(N_mesh*2+1))
plt.figure(f)
flatLinePlot(X, scalar_flux)
plt.title('Infinte Med')
plt.xlabel('Distance')
plt.ylabel('Scalar Flux')
plt.show()

f+=1
plt.figure(f)
plt.title('Infinte Med')
plt.xlabel('Distance')
plt.ylabel('Current')
plt.ylim([-1,1])
flatLinePlot(X, current)
plt.show()
#'''

