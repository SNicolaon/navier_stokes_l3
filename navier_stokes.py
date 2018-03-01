# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 17:25:58 2018

@author: Sacha
"""

from matplotlib import pyplot as plt
import math
import numpy as np

#parametre:
t=.5 #temps
lx=10 #largeur
ly=10 #hauteur

g=9.81 #pas utilisé , a voir

patm=101325
rho=1.225 # air : 1.225
nu = 15.6* 10**(-6) #visco cinema air : 15.6* 10**(-6)
F=300 #difference de pression tq v moyen=150m/s : 300


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

p = np.ones([ni,nj]) #pression #CI : P=1
p=patm*p



#fonctions:



def calcul_p(u,v,p): #p(t) tq div( v(t) ) =0 
    pn = np.zeros_like(p)
    
    pn=p.copy()
    
    pn[1:-1, 1:-1] = (((p[1:-1, 2:] + p[1:-1, 0:-2]) * dy**2 + (p[2:, 1:-1] +
                      p[0:-2, 1:-1]) * dx**2) / (2 * (dx**2 + dy**2)) - 
                          dx**2 * dy**2 / (2 * (dx**2 + dy**2)) *
                          (  (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
                        (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -  
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 * ((u[2:, 1:-1] -
                    u[0:-2, 1:-1]) / (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))- 
                          ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2)) ))
    
    
    pn[1:-1, -1] = (((p[1:-1, 0] + p[1:-1, -2])* dy**2 +
                        (p[2:, -1] + p[0:-2, -1]) * dx**2) /
                       (2 * (dx**2 + dy**2)) -
                       dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * ( (rho * (1 / dt *
                                       ((u[1:-1, 0] - u[1:-1,-2]) / (2 * dx) +
                                    (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) -
                          ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 -
                          2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) *
                               (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) -
                          ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2))   ))
    
    
    pn[1:-1, 0] = (((p[1:-1, 1] + p[1:-1, -1])* dy**2 +
                       (p[2:, 0] + p[0:-2, 0]) * dx**2) /
                      (2 * (dx**2 + dy**2)) -
                      dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * ( (rho * (1 / dt *
                                      ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) +
                                   (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) -
                         ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx))**2 -
                         2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) *
                              (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx))-
                         ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))**2))    ))
    #CL:
    
    pn[ni-1,:]=pn[ni-2,:]  # dp/dy =0 et p=0 sur les murs
    pn[0,:]=pn[1,:]    #aux murs : remarque : on a pas pn[1:-1, 0] et pn[1:-1, nj-1] : a calculer separement
    
    
    return(pn)



def calcul_v(u,v,p): #calcul u,v (t+1)
    un= np.zeros_like(u)
    un=u.copy()
    vn= np.zeros_like(v)
    vn=v.copy()
    
    
    un[1:-1, 1:-1] = (u[1:-1, 1:-1]- u[1:-1, 1:-1] * dt / dx *
                      (u[1:-1, 1:-1] - u[1:-1, 0:-2]) - v[1:-1, 1:-1] * dt / dy *
                      (u[1:-1, 1:-1] - u[0:-2, 1:-1]) - 
                         dt / (2 * rho * dx) * (p[1:-1, 2:] -
                     p[1:-1, 0:-2]) + nu *
                    (dt / dx**2 * (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, 0:-2]) + 
                         dt / dy**2 * (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1]))) + F*dt
    
    vn[1:-1,1:-1] = (v[1:-1, 1:-1] - u[1:-1, 1:-1] * dt / dx * (v[1:-1, 1:-1] - v[1:-1, 0:-2]) -
                      v[1:-1, 1:-1] * dt / dy *(v[1:-1, 1:-1] - v[0:-2, 1:-1]) - 
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu *
                        (dt / dx**2 * (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, 0:-2]) + 
                        dt / dy**2 * (v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[0:-2, 1:-1])))
    
    
    un[1:-1, -1] = (u[1:-1, -1] - u[1:-1, -1] * dt / dx * 
                  (u[1:-1, -1] - u[1:-1, -2]) -
                   v[1:-1, -1] * dt / dy * 
                  (u[1:-1, -1] - u[0:-2, -1]) -
                   dt / (2 * rho * dx) *
                  (p[1:-1, 0] - p[1:-1, -2]) + 
                   nu * (dt / dx**2 * 
                  (u[1:-1, 0] - 2 * u[1:-1,-1] + u[1:-1, -2]) +
                   dt / dy**2 * 
                  (u[2:, -1] - 2 * u[1:-1, -1] + u[0:-2, -1])) + F * dt)
    
    un[1:-1, 0] = (u[1:-1, 0] - u[1:-1, 0] * dt / dx *
                 (u[1:-1, 0] - u[1:-1, -1]) -
                  v[1:-1, 0] * dt / dy * 
                 (u[1:-1, 0] - u[0:-2, 0]) - 
                  dt / (2 * rho * dx) * 
                 (p[1:-1, 1] - p[1:-1, -1]) + 
                  nu * (dt / dx**2 * 
                 (u[1:-1, 1] - 2 * u[1:-1, 0] + u[1:-1, -1]) +
                  dt / dy**2 *
                 (u[2:, 0] - 2 * u[1:-1, 0] + u[0:-2, 0])) + F * dt)
    
    vn[1:-1, -1] = (v[1:-1, -1] - u[1:-1, -1] * dt / dx *
                  (v[1:-1, -1] - v[1:-1, -2]) - 
                   v[1:-1, -1] * dt / dy *
                  (v[1:-1, -1] - v[0:-2, -1]) -
                   dt / (2 * rho * dy) * 
                  (p[2:, -1] - p[0:-2, -1]) +
                   nu * (dt / dx**2 *
                  (v[1:-1, 0] - 2 * v[1:-1, -1] + v[1:-1, -2]) +
                   dt / dy**2 *
                  (v[2:, -1] - 2 * v[1:-1, -1] + v[0:-2, -1])))
    
    
    vn[1:-1, 0] = (v[1:-1, 0] - u[1:-1, 0] * dt / dx *
                 (v[1:-1, 0] - v[1:-1, -1]) -
                  v[1:-1, 0] * dt / dy *
                 (v[1:-1, 0] - v[0:-2, 0]) -
                  dt / (2 * rho * dy) * 
                 (p[2:, 0] - p[0:-2, 0]) +
                  nu * (dt / dx**2 * 
                 (v[1:-1, 1] - 2 * v[1:-1, 0] + v[1:-1, -1]) +
                  dt / dy**2 * 
                 (v[2:, 0] - 2 * v[1:-1, 0] + v[0:-2, 0])))
    
    
    
    
    #CL:
    un[ni-1,:]=0.
    vn[0,:]=0.
    vn[ni-1,:]=0. # dv/dy =0 et v=0 sur les murs
    un[0,:]=0. #aux murs : remarque : on a pas v/u n[1:-1, 0] et v/u n[1:-1, nj-1] : a calculer separement
    

    return(un,vn)





#calcul:
for ti in range(it):
    for k in range(lapla):
        p=calcul_p(u,v,p)
    
    u,v=calcul_v(u,v,p)
    



#affichage:
fig = plt.figure()

x = np.linspace(0, lx, nj)
y = np.linspace(0, ly, ni)
X, Y = np.meshgrid(x, y)

plt.quiver(X,Y, u, v) 

plt.contourf(X, Y, p, alpha=0.7)  
plt.colorbar()
