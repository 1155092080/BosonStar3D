#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import time
import gc
from scipy.interpolate import interp1d


# In[2]:

#Grid Variables
N = 300
half = 150
a = -128.0
b = 128.0

dt = 10.0
tMax = 1000.1
outputEvery = 100.0;

#Physical variables
epsilon = 1.0
k1 = 4.32747
G = 0.344369
M = 6.53504
Ra = 44614.3

#Working area
rstatic, rhostatic = np.loadtxt("stable_star.txt", unpack=True)
sourceFile = open("demo.txt", 'w')

rstatic /= 6.77193e-6 * 100
rstatic /= Ra
rhostatic /= 1.619e-18
rhostatic *= 1000
rhostatic /= (M*1.9891e30)
psistatic = np.sqrt(rhostatic)
psistatic *= np.power(Ra, 3/2)
    
#Init
x_ = np.linspace(a,b,N)
y_ = np.linspace(a,b,N)
z_ = np.linspace(a,b,N)
x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

x2 = x*x
y2 = y*y
z2 = z*z
r2 = (x2+y2+z2)

muLx = np.linspace(-half, half-1, N)
muLx *= 2.0*np.pi/(b-a)
muLy = np.linspace(-half, half-1, N)
muLy *= 2.0*np.pi/(b-a)
muLz = np.linspace(-half, half-1, N)
muLz *= 2.0*np.pi/(b-a)
muLxmesh, muLymesh, muLzmesh = np.meshgrid(muLx, muLy, muLz, indexing='ij')

muLr2 = muLxmesh*muLxmesh+muLymesh*muLymesh+muLzmesh*muLzmesh

psi = np.interp(np.sqrt(r2), rstatic, psistatic, left=psistatic[0], right=0.0)
psi = psi + 0.0*1j

phi3D = np.empty((N,N,N))
    
count = 0
t = 0.0
nextOutput = 0.0;

# In[3]:


k1dtdivtwoepsilon = k1*dt/2.0/epsilon

del x2
del y2
del z2
del x
del y
del z
del muLxmesh, muLymesh, muLzmesh
gc.collect()

def findPotential(rdata):
    global phi3D
    tolerance = 3e-8
    xdata = np.sqrt(r2[half][half][half:])
    N2 = len(rdata)
    dx = np.diff(xdata, prepend=xdata[0])

    temp = np.empty(N2)
    temp[0] = (4.0/3.0)*np.pi*xdata[0]*xdata[0]*xdata[0]*rdata[0]
    temp[1] = 9.0*temp[0] + (4.0/3.0)*np.pi*dx[1]*xdata[1]*xdata[1]*rdata[1];
    for i in range(2,N2):
        temp[i] = temp[i-2] + 4.0 * np.pi * (dx[i]/3.0) * \
                     (rdata[i-2] * xdata[i-2] * xdata[i-2] + \
                     4.0*rdata[i-1]*xdata[i-1]*xdata[i-1] + \
                     rdata[i] * xdata[i]*xdata[i])

    phip = temp / (xdata*xdata)

    phi = np.empty(N2)
    phi[N2-1] = 0.0;

    for i in range(N2-2,-1,-1):
        phi[i] = phi[i+1]-(phip[i]+phip[i+1])*dx[i]/2.0

    #Relax
    phiNew = np.empty(N2)
    error = np.empty(N2)

    while(True):

        phiNew[0] = 0.5 * (phi[1] + phi[1]) - 2.0*np.pi*dx[0]*dx[0]*rdata[0]
        for j in range(1,N2-1):
            phiNew[j] = 0.5 * (phi[j+1] + phi[j-1]) + 0.5 * dx[j]/xdata[j] * (phi[j+1] - phi[j-1]) - 2.0*np.pi*dx[j]*dx[j]*rdata[j]

        phiNew[N2-1] = phi[N2-1]

        for j in range(N2-1):
            error[j] = (phiNew[j] - phi[j])/phi[j]

        error[N2-1] = 0.0

        phi = phiNew

        if(error.any() > tolerance):
            continue
        else:
            #Get V(x,y,z) here
            phiInt = interp1d(xdata, phi*G, kind='cubic', bounds_error=False, fill_value=(phi[0]*G,0.0))
            phi3D = phiInt(np.sqrt(r2))
            break

    return


# In[4]:


while(True):

    if(t >= nextOutput):
        plt.figure()
        plt.plot(x_, np.real(psi[half][half][:]*np.conj(psi[half][half][:])))
        plt.ylim(-0.00002,0.0002)
        plt.savefig("ttt"+str(count)+".pdf")
        #plt.show()
        nextOutput += outputEvery;
        count += 1
    
    start = time.time()
    findPotential(np.real(psi[half][half][half:]*np.conj(psi[half][half][half:])))

    tmp = phi3D*dt/2.0/epsilon + k1dtdivtwoepsilon*psi*np.conj(psi)

    psiStar = (np.cos(tmp) - np.sin(tmp)*1j)*psi
    end1 = time.time()
    
    print("Action 1:", end1-start)

    psiFourier = np.fft.fftn(psiStar, axes=(-3,-2,-1))
    psiFourier = np.fft.fftshift(psiFourier)
    end2 = time.time()
    
    print("Action 2:", end2-start)
    
    psiFourier = (np.cos(epsilon*dt*muLr2/2.0) - np.sin(epsilon*dt*muLr2/2.0)*1j)*psiFourier
    end3 = time.time()
    
    print("Action 3:", end3-start)
    
    psiFourier = np.fft.ifftshift(psiFourier)
    psiStarStar = np.fft.ifftn(psiFourier, axes=(-3,-2,-1))
    end4 = time.time()
    
    print("Action 4:", end4-start)

    findPotential(np.real(psiStarStar[half][half][half:]*np.conj(psiStarStar[half][half][half:])))
    
    tmp = phi3D*dt/2.0/epsilon + k1dtdivtwoepsilon*psiStarStar*np.conj(psiStarStar)
    psi = (np.cos(tmp) - np.sin(tmp)*1j)*psiStarStar
    
    end = time.time()
    print(t, end-start)
    
    t += dt
    if(t >= tMax):
        break

    print(np.real(psi[half][half][half]*np.conj(psi[half][half][half])), file = sourceFile)
    
sourceFile.close()

# In[ ]:




