import numpy as np
import matplotlib.pyplot as plt

def norm(vec):
    for i in range(vec.shape[1]):
        vec[:,i] = vec[:,i] / max(vec[:,i] )
    
    return(vec)

dt = np.array([.1, 0.075, 0.05, 0.025, 0.0125, 0.01, 0.0075, 0.005, 0.001])
N_dt = 5

#np.savez('output_dt_{0}'.format(dt[i]), sfEuler=sfEuler, sfMB=sfMB, specRadsMB=specRadsMB, specRadsEuler=specRadsEuler)
errorMB = np.zeros(N_dt)
errorEuler = np.zeros(N_dt)

fig, axs = plt.subplots(N_dt)
for i in range(N_dt):
    x = np.linspace(0, 1, int(200))

    print()
    print()
    print('>>>>>>>>>>>>>>> dt: {0} <<<<<<<<<<<<<<<<'.format(dt[i]) )
    print()
    with np.load('detData/output_dt_{0}.npz'.format(dt[i])) as detData:
        sfEuler = detData['sfEuler']
        sfMB = detData['sfMB']
    with np.load('mcdcData/output_mcdc_dt_{0}.npz'.format(dt[i])) as mcData:
        sfRef = mcData['sfRef']

    mid = int(sfRef.shape[1]/2)
    print()
    print('shape MB:    {0}'.format(sfMB.shape))
    print('shape Euler: {0}'.format(sfEuler.shape))
    print('shape Ref:   {0}'.format(sfRef.shape))
    print()
    print('midval MB    {0}'.format(sfMB[mid, -1]))
    print('midval E     {0}'.format(sfEuler[mid, -1]))
    print('midval Ref   {0}'.format(sfRef[mid, -1]))
    print()

    sfEuler = norm(sfEuler)
    sfMB    = norm(sfMB)
    sfRef   = norm(sfRef)

    print()
    print('>>>after norming')
    print('midval MB    {0}'.format(sfMB[mid, -1]))
    print('midval E     {0}'.format(sfEuler[mid, -1]))
    print('midval Ref   {0}'.format(sfRef[mid, -1]))
    print()

    axs[i].plot(x, sfMB[:, -2])
    axs[i].plot(x, sfEuler[:, -2])
    axs[i].plot(x, sfRef[:, -1])
    axs[i].set_title('dt: {0}'.format(dt[i]))
    for ax in axs.flat:
        ax.label_outer()

    errorMB[i] =    np.linalg.norm((sfMB[:,-1]-sfRef[:,-1]) / sfRef[:,-1])
    errorEuler[i] = np.linalg.norm((sfEuler[:,-1]-sfRef[:,-1]) / sfRef[:,-1])


plt.show()

print()
print(errorMB)
print()
print(errorEuler)
print()

plt.figure()
plt.loglog(dt[:N_dt], errorMB, label='MB')
plt.loglog(dt[:N_dt], errorEuler, label='E')
plt.legend()
plt.show()
