# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 17:25:58 2018

@author: Sacha
"""

from matplotlib import pyplot as plt
#import math
import numpy as np
from scipy import misc

###parametre:
t=2**(-2) #temps
lx=10 #largeur
ly=10 #hauteur

#g=9.81 #pas utilisé , a voir

patm=105 #101325                                              #cause overflow
rho=1.225 # air : 1.225
nu = 0.84  #visco cinema air : 15.6* 10**(-6) ou  2**(-16)      #cause overflow
#F=1 #force appliqué au fluide(difference de pression)      #cause overflow 
vm=5 #vitesse dentree et de sorti

#discretisation :
nj=101 #largeur                                              #overflow en 200*100
ni=51 #hauteur
dt=2**(-10)  #pas de temps =0.001 mais pas de 1/10 (non exacte => overflow plus fcilement)
lapla=20 #iteration du poisson (calcul pression)



#intervalles:
dx= lx / (nj-1)
dy= ly / (ni-1)
it= int((t / dt ))+1  #nombre d'iteration temps


#variable:
u = np.ones([ni,nj],dtype=np.float64) #vitesse sur x
v = np.zeros([ni,nj],dtype=np.float64) # sur y

u=vm*u
p = np.ones([ni,nj],dtype=np.float64) #pression #CI : P=atm
p=patm*p



###fonctions:



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
    
    """
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
    """
    
    #CL:
    
    pn[ni-1,:]=pn[ni-2,:]  # dp/dy =0 et p=0 sur les murs
    pn[0,:]=pn[1,:]    #aux murs : remarque : on a pas pn[1:-1, 0] et pn[1:-1, nj-1] : a calculer separement
    
    pn=pn*Mtot
    pn=pn  + 0*Mint + Mfront*(    np.concatenate( (np.zeros([ni,1]) , np.concatenate( (np.zeros([1,nj-2]) , (pn[1:-1, 0:-2]+pn[1:-1, 2:] +pn[ 0:-2 , 1:-1]+pn[ 2: , 1:-1])/2 ,  np.zeros([1,nj-2])), axis=0 ),  np.zeros([ni,1])), axis=1 )     )
    
    pn[1:-1, 0]=patm
    pn[1:-1, -1]=patm
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
                         dt / dy**2 * (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1]))) #+ F*dt 
    
    vn[1:-1,1:-1] = (v[1:-1, 1:-1] - u[1:-1, 1:-1] * dt / dx * (v[1:-1, 1:-1] - v[1:-1, 0:-2]) -
                      v[1:-1, 1:-1] * dt / dy *(v[1:-1, 1:-1] - v[0:-2, 1:-1]) - 
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu *
                        (dt / dx**2 * (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, 0:-2]) + 
                        dt / dy**2 * (v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[0:-2, 1:-1])))
    
    """
    un[1:-1, -1] = (u[1:-1, -1] - u[1:-1, -1] * dt / dx * 
                  (u[1:-1, -1] - u[1:-1, -2]) -
                   v[1:-1, -1] * dt / dy * 
                  (u[1:-1, -1] - u[0:-2, -1]) -
                   dt / (2 * rho * dx) *
                  (p[1:-1, 0] - p[1:-1, -2]) + 
                   nu * (dt / dx**2 * 
                  (u[1:-1, 0] - 2 * u[1:-1,-1] + u[1:-1, -2]) +
                   dt / dy**2 * 
                  (u[2:, -1] - 2 * u[1:-1, -1] + u[0:-2, -1])) )#+ F * dt)
    
    un[1:-1, 0] = (u[1:-1, 0] - u[1:-1, 0] * dt / dx *
                 (u[1:-1, 0] - u[1:-1, -1]) -
                  v[1:-1, 0] * dt / dy * 
                 (u[1:-1, 0] - u[0:-2, 0]) - 
                  dt / (2 * rho * dx) * 
                 (p[1:-1, 1] - p[1:-1, -1]) + 
                  nu * (dt / dx**2 * 
                 (u[1:-1, 1] - 2 * u[1:-1, 0] + u[1:-1, -1]) +
                  dt / dy**2 *
                 (u[2:, 0] - 2 * u[1:-1, 0] + u[0:-2, 0])) )#+ F * dt)
    
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
    
    """
    
    
    #CL:
    un[ni-1,:]=0.
    vn[0,:]=0.
    vn[ni-1,:]=0. # dv/dy =0 et v=0 sur les murs
    un[0,:]=0. #aux murs : remarque : on a pas v/u n[1:-1, 0] et v/u n[1:-1, nj-1] : a calculer separement
    
    un[1:-1,0]=vm
    un[1:-1,-1]=vm
    vn[1:-1,0]=0
    vn[1:-1,-1]=0
    
    vn=vn*Mtot
    un=un*Mtot
    
    return(un,vn)


###Condition limites

tab = misc.imread('wing.png') # 41*41*4 array   #il semblerait que j'ai codé la verticalité a l'envers => symetrie verticale sur l'image
#[y,x,couleur]
#white : (255)(255)(255)   (255)
#black : ( 0 )( 0 )( 0 )   (255)
#red:    (255)( 0 )( 0 )   (255)

Mfront=np.zeros([ni,nj],np.float64) #black = 1  white = 0 red=0
Mtot=np.ones([ni,nj],np.float64)#black = 0  white = 1 red=0
Mint=np.zeros([ni,nj],np.float64)#black = 0  white = 0 red=1 #innutile pour le moment
#Mtot est le profil entier , Mint que l'interieur (sans frontiere)

for i in range(ni):
    for j in range(nj):
        if(tab[i,j,0]==0): #black
            Mfront[i,j]=1
            Mtot[i,j]=0
            Mint[i,j]=0
        if(tab[i,j,0]==255):
            if(tab[i,j,1]==0): #red
                Mfront[i,j]=0
                Mtot[i,j]=0
                Mint[i,j]=1




###calcul:
for ti in range(it):
    for k in range(lapla):
        p=calcul_p(u,v,p)
    
    u,v=calcul_v(u,v,p)
    



###affichage:
fig = plt.figure()

x = np.linspace(0, lx, nj)
y = np.linspace(0, ly, ni)
X, Y = np.meshgrid(x, y)

plt.quiver(X,Y, u, v) 

plt.contourf(X, Y, p, alpha=0.7)  
plt.colorbar()
