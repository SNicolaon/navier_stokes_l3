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
t=10**(-2) #temps
lx=10 #largeur
ly=10 #hauteur


patm=101325   #101325                                           
rho=1.225 # air : 1.225
nu = 1  #visco cinema air : 15.6* 10**(-6) ou  2**(-16)      
vm=30 #vitesse dentree et de sorti

#discretisation :
nj=51 #largeur                                            
ni=51 #hauteur
dt=10**(-5)  #pas de temps =0.001 la base 
lapla=10 #iteration du poisson (calcul pression)



#intervalles:
dx= lx / (nj-1)
dy= ly / (ni-1)
it= int((t / dt ))+1  #nombre d'iteration temps


#variable:
u = np.ones([ni,nj],dtype=np.float64) #vitesse sur x
v = np.zeros([ni,nj],dtype=np.float64) # sur y
p = np.ones([ni,nj],dtype=np.float64) #pression 

u=vm*u

u[:,0]=vm #soit ca soit celui qu'es dans la boucle V , lui converge mieux mais faut troouver un moyen de le faire parabolique sinon trou de pression
u[:,-1]=vm

u[0,:]=0 
u[ni-1,:]=0 #vm ou 0 si y a un mur (conduite)

p=p*patm 


###fonctions:


CFL=vm*dt/dx
print(CFL)

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
    
  
    
    #CL:
    """
    pn[ni-1,:]=pn[ni-2,:]  # dp/dy =0 sur les murs
    pn[0,:]= pn[1,:]   
    
    
    
    pn[1:-1, -1]=pn[1:-1, -2] #pas une bonne CL mais honntement pour ce que ca change
    """
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
                     p[1:-1, 0:-2]) +
                        nu *  (dt / dx**2 * (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, 0:-2]) + 
                         dt / dy**2 * (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1]))) 
    
    vn[1:-1,1:-1] = (v[1:-1, 1:-1] - u[1:-1, 1:-1] * dt / dx * (v[1:-1, 1:-1] - v[1:-1, 0:-2]) -
                      v[1:-1, 1:-1] * dt / dy *(v[1:-1, 1:-1] - v[0:-2, 1:-1]) - 
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) + 
                        nu * (dt / dx**2 * (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, 0:-2]) + 
                        dt / dy**2 * (v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[0:-2, 1:-1])))
    

    
    
    #CL:

    vn=vn*Mtot #objet
    un=un*Mtot
    
    
    """
    un[:,0]=un[:,1] #dV/dx =0 a l'entree et sortie ; pas d'acceleration
    un[:,-1]=un[:,-2]
    """
    return(un,vn)


###Condition limites

tab = misc.imread('car.png') # 41*41*4 array  
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
    print(100*ti/it,'%')
    



###affichage:
fig = plt.figure()

x = np.linspace(0, lx, nj)
y = np.linspace(0, ly, ni)
X, Y = np.meshgrid(x, y)

plt.quiver(X,Y, u, v) 

plt.contourf(X, Y, p, alpha=0.7)  
plt.colorbar()

plt.show()
