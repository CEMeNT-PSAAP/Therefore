# method of manufactured soultions plotting

import numpy as np
import matplotlib.pyplot as plt


a = 7.9
b = 2
c = 4
xsec = 0.25
xsecScatter = 0.1
v = 2
Nx = 30

def Q(x,t,mu,l):
    Qm = 2*b*np.exp(l*t)*(1/v*x*t + mu + xsec*x) + 2*mu**2 + 2*xsec*c*x*mu + xsec*a - xsecScatter*c/3*x 
    return(Qm)

def af(x,t,mu,l):
    afm = a+ b*x*np.exp(l*t) + c*mu*x

x = np.linspace(0, 100, Nx)
print(x)

plt.figure()

afP = []
QP = []
for i in range(Nx):
    afP.append(af(x[i],0.5,0.717,-0.12))
    QP.append(Q(x[i],0.5,0.717,-0.12))



plt.plot(x, afP)
plt.plot(x, QP)
plt.show()
