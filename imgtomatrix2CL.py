# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 10:13:29 2018

@author: Sacha
"""


#white : (255)(255)(255)   (255)
#black : ( 0 )( 0 )( 0 )   (255)
#red:    (255)( 0 )( 0 )   (255)


import numpy as np
from scipy import misc



tab = misc.imread('profil.png') # 41*41*4 array
#[y,x,couleur]
#white : (255)(255)(255)   (255)
#black : ( 0 )( 0 )( 0 )   (255)
#red:    (255)( 0 )( 0 )   (255)

Mfront=np.zeros([41,41]) #black = 1  white = 0 red=0
Mtot=np.ones([41,41])#black = 0  white = 1 red=0
Mint=np.ones([41,41])#black = 1  white = 1 red=0 #innutile pour le moment
#Mtot est le profil entier , Mint que l'interieur (sans frontiere)

for i in range(41):
    for j in range(41):
        if(tab[i,j,0]==0):
            print('black')
            Mfront[i,j]=1
            Mtot[i,j]=0
            Mint[i,j]=1
        if(tab[i,j,0]==255):
            if(tab[i,j,1]==0):
                print('red')
                Mfront[i,j]=0
                Mtot[i,j]=0
                Mint[i,j]=0




#calcul

V=np.ones([41,41])
V=V*Mtot #pour Vx et Vy




P=np.ones([41,41])
P=P*Mtot

Ptemp1=(P[1:-1, 0:-2]+P[1:-1, 2:] +P[ 0:-2 , 1:-1]+P[ 2: , 1:-1])/2  #taille 39*39
b=np.zeros([1,39])
Ptemp2=np.concatenate( (b , Ptemp1,  b), axis=0 ) #39*41
c=np.zeros([41,1])
Ptemp3=np.concatenate( (c , Ptemp2,  c), axis=1 ) #41*41

P=P+ Mfront*Ptemp3
