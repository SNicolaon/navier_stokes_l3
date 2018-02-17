# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 17:25:58 2018

@author: Sacha
"""

from matplotlib import pyplot as plt
import math
import numpy as np
import imageio

#parametre:
t=2.0 #temps
lx=4 #largeur
ly=4 #hauteur

rho=1
g=9.81
nu = 0.1

#discretisation :
nj=41 #largeur
ni=41 #hauteur
dt=0.001  #pas de temps
lapla=50 #iteration du laplacien (calcul pression)

#intervalles:
dx= lx / (nj-1)
dy= ly / (ni-1)
it= int((t / dt ))+1  #nombre d'iteration temps

#variable:
u = np.zeros([ni,nj]) #vitesse sur x
v = np.zeros([ni,nj]) # sur y
p = np.zeros([ni,nj]) #pression





#fonctions:

"""
def cond_lim(u,v,p,im):
    if (im[i,j]== blabla) :
        u[i,j]= 0
        v[i,j]= 0



"""


def calcul_p(u,v,p): #p(t) tq div( v(t) ) =0 
    pn = np.zeros_like(p)
    
    pn=p.copy()
    
    pn[1:-1, 1:-1] = (((p[1:-1, 2:] + p[1:-1, 0:-2]) * dy**2 + (p[2:, 1:-1] + p[0:-2, 1:-1]) * dx**2) / (2 * (dx**2 + dy**2)) - \
                          dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * (  (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -  \
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))- \
                          ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2)) ))
    
    #CL:
    
    pn[ni-1,:]=0.
    pn[0,:]=pn[1,:]
    pn[:,0]=pn[:,1]
    pn[:,nj-1]=pn[:,nj-2]
    
    return(pn)



def calcul_v(u,v,p): #calcul u,v (t+1)
    un= np.zeros_like(u)
    un=u.copy()
    vn= np.zeros_like(v)
    vn=v.copy()
    
    
    un[1:-1, 1:-1] = (u[1:-1, 1:-1]- u[1:-1, 1:-1] * dt / dx * (u[1:-1, 1:-1] - u[1:-1, 0:-2]) - v[1:-1, 1:-1] * dt / dy *(u[1:-1, 1:-1] - u[0:-2, 1:-1]) - \
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) + nu * (dt / dx**2 * (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, 0:-2]) + \
                         dt / dy**2 * (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1])))
    
    vn[1:-1,1:-1] = (v[1:-1, 1:-1] - u[1:-1, 1:-1] * dt / dx * (v[1:-1, 1:-1] - v[1:-1, 0:-2]) - v[1:-1, 1:-1] * dt / dy *(v[1:-1, 1:-1] - v[0:-2, 1:-1]) - \
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu * (dt / dx**2 * (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, 0:-2]) + \
                        dt / dy**2 * (v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[0:-2, 1:-1])))
    
    #CL:
    un[ni-1,:]=1.
    vn[0,:]=0.
    
    un[0,:]=0.
    vn[ni-1,:]=0.
    
    un[:,0]=0.
    vn[:,0]=0.
    
    un[:,nj-1]=0.
    vn[:,nj-1]=0.
    
    un[ni-1,:]=1.
    return(un,vn)


"""
#conditions limites:
im = imageio.imread('profil.png') #importe l'image en forme de tableau
print(im.shape)
print(im)

cond_lim()
"""





#calcul:
for ti in range(it):
    for k in range(lapla):
        p=calcul_p(u,v,p)
    
    u,v=calcul_v(u,v,p)
    


print(1)

#affichage:
fig = plt.figure()

x = np.linspace(0, lx, nj)
y = np.linspace(0, ly, ni)
X, Y = np.meshgrid(x, y)

plt.quiver(X,Y, u, v) 

plt.contourf(X, Y, p, alpha=0.7)  
plt.colorbar()

