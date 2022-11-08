#!/usr/bin/env python
# coding: utf-8

# In[1]:

import sys
import numpy as np
import matplotlib

matplotlib.use('Agg')

from matplotlib import pyplot as plt
import time
import gc
from scipy.interpolate import interp1d

# In[2]:

#Grid Variables
N = 800
a = -10.0
b = 10.0

#Need to add this part (Read from argument)
dt = float(sys.argv[1])

#Physical variables
epsilon = 1.0
k1 = 563.835
G = 44.8685
M = 0.47386
Ra = 390012.0
omega = 4.8614

#Working area
rstatic, rhostatic, phistatic = np.loadtxt("./Temp/tmp001.txt", unpack=True)
sourceFile = open("./Temp/tmp002.txt", 'w')

rstatic /= 6.77193e-6 * 100
rstatic /= Ra
rhostatic /= 1.619e-18
rhostatic *= 1000
rhostatic /= (M*1.9891e30)
psistatic = np.sqrt(rhostatic)
psistatic *= np.power(Ra, 3/2)
phistatic *= (1e-7/5.59429e-55)
phistatic *= (5.02788e-34*1000)
phistatic *= np.power(Ra, -2)
phistatic *= np.power(omega, -2)
dt /= 2.03017e5
dt *= omega

#Init
x_ = np.linspace(a,b,N)
y_ = np.linspace(a,b,N)
z_ = np.linspace(a,b,N)
x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

x2 = x*x
y2 = y*y
z2 = z*z
r2 = (x2+y2+z2)

muLx = np.linspace(-N/2, N/2-1, N)
muLx *= 2.0*np.pi/(b-a)
muLy = np.linspace(-N/2, N/2-1, N)
muLy *= 2.0*np.pi/(b-a)
muLz = np.linspace(-N/2, N/2-1, N)
muLz *= 2.0*np.pi/(b-a)
muLxmesh, muLymesh, muLzmesh = np.meshgrid(muLx, muLy, muLz, indexing='ij')

muLr2 = muLxmesh*muLxmesh+muLymesh*muLymesh+muLzmesh*muLzmesh

psi = np.interp(np.sqrt(r2), rstatic, psistatic, left=psistatic[0], right=0.0)
#rhoInt = interp1d(rstatic, psistatic, kind='cubic', bounds_error=False, fill_value=(psistatic[0],0.0))
#psi = rhoInt(np.sqrt(r2))
psi = psi + 0.0*1j

#phi3D = np.empty((N,N,N))
#phiInt = interp1d(rstatic, phistatic, kind='cubic', bounds_error=False, fill_value=(phistatic[0],0.0))
#phi3D = phiInt(np.sqrt(r2))
phi3D = np.interp(np.sqrt(r2), rstatic, phistatic, left=phistatic[0], right=0.0)


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

# In[4]:
    
plt.figure()
plt.plot(x_, np.real(psi[400][400][:]*np.conj(psi[400][400][:])))
#plt.ylim(-0.00002,0.0002)
plt.savefig("test.pdf")

start = time.time()
#findPotential(np.real(psi[400][400][400:]*np.conj(psi[400][400][400:])))

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

#findPotential(np.real(psiStarStar[400][400][400:]*np.conj(psiStarStar[400][400][400:])))

tmp = phi3D*dt/2.0/epsilon + k1dtdivtwoepsilon*psiStarStar*np.conj(psiStarStar)
psi = (np.cos(tmp) - np.sin(tmp)*1j)*psiStarStar

end = time.time()
print(end-start)

DMdensity = np.interp(rstatic, np.sqrt(r2[400][400][400:]), np.real(psi[400][400][400:]*np.conj(psi[400][400][400:])), left=np.real(psi[400][400][400]*np.conj(psi[400][400][400])), right=0.0)

rstatic *= 6.77193e-6 * 100
rstatic *= Ra
DMdensity /= (Ra*Ra*Ra)
DMdensity *= 1.619e-18
DMdensity /= 1000
DMdensity *= (M*1.9891e30)

for i in range(len(rstatic)):
    print(str(rstatic[i])+"\t"+str(DMdensity[i]), file = sourceFile)
    
sourceFile.close()

# In[ ]: