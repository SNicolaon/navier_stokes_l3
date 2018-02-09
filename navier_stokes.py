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
t=10 #temps
lx=10 #largeur
ly=10 #hauteur

rho=1
g=9.81
nu = 0.1

#discretisation :
nj=50 #largeur
ni=50 #hauteur
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

f = np.zeros([ni,nj]) #calcul de poisson

#affichage:
U = np.zeros([ni,nj,it]) #les U a tout les temps
V = np.zeros([ni,nj,it]) # les V
P = np.zeros([ni,nj,it]) # les P

#fonctions:

"""
def con_lim(u,v,p,im):
    if im[i,j]:
        u[i,j]= 0
        v[i,j]= 0



"""

def calcul_p(u,v,p): #p(t) tq div( v(t) ) =0 
    pn = np.empty_like(p)
    for j in range(1,nj-1): #on connait les CL des 2cotes 
        for i in range(1,ni-1):
            pn[i,j]= ((((p[i+1,j] + p[i-1,j])*dy**2 + (p[i,j+1] + p[i,j-1])*dx**2) / (2*(dx**2 + dy**2)) ) - ((rho* dx**2 * dy**2 )/ 2(dx**2 + dy**2))*( (((u[i+1,j]-u[i-1,j])/(2*dx)) + ((v[i,j+1]-v[i,j-1])/(2*dy)))/dt -  ((u[i+1,j]-u[i-1,j])/(2*dx))**2 - ((v[i,j+1]-v[i,j-1])/(2*dy))**2  - (1/(2*dx*dy))*(u[i,j+1]-u[i,j-1])*(v[i+1,j]-v[i-1,j])))
    p = pn.copy()
    return(p)



def calcul_v(u,v,p): #calcul u,v (t+1)
    un= np.empty_like(u)
    vn= np.empty_like(v)
    for j in range(1,nj-1):
        for i in range(1,ni-1):
            un[i,j]=u[i,j]- v[i,j]*u[i,j]*(dt**2 /(dy*dx))* (u[i,j] + u[i-1,j])*(u[i,j] + u[i,j-1]) - (dt / (2*rho*dx)) *(p[i+1,j] - p[i-1,j]) + nu*dt*(((u[i+1,j] -2* u[i,j] + u[i-1,j]  )/dx**2)  + ( ( u[i,j+1] -2* u[i,j] + u[i,j-1]     )/dy**2))
            vn[i,j]=v[i,j]- v[i,j]*u[i,j]*(dt**2 /(dy*dx))* (v[i,j] + v[i-1,j])*(v[i,j] + v[i,j-1]) - (dt / (2*rho*dy)) *(p[i,j+1] - p[i,j-1]) + nu*dt*(((v[i+1,j] -2* v[i,j] + v[i-1,j]  )/dx**2)  + ( ( v[i,j+1] -2* v[i,j] + v[i,j-1]     )/dy**2))
    u = un.copy()
    v = vn.copy()
    return(u,v)

"""
#conditions limites:
im = imageio.imread('profil.png') #importe l'image en forme de tableau
print(im.shape)
print(im)
"""

"""
#grid:
x = np.linspace(0, lx, nj)
y = np.linspace(0, ly, ni)
X, Y = np.meshgrid(x, y)

"""




"""
u[1:] - u[0:-1]

{chaque case de u a partir de 1} - {chaque case de u de 0 a taille(u)-1 }

optimise le temps de calcul (genre de 1s pour un calcul simple)
"""




#calcul:
for ti in range(it):
    for k in range(lapla):
        p=calcul_p(u,v,p)
    
    u,v=calcul_v(u,v,p)
    
    U[:,:,ti]=u[:,:]
    V[:,:,ti]=v[:,:]
    P[:,:,ti]=p[:,:]







