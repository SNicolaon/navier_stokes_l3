# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 10:13:29 2018

@author: Sacha
"""


#white : (255)(255)(255)   (255)
#black : ( 0 )( 0 )( 0 )   (255)
#red:    (255)( 0 )( 0 )   (255)

#transforme ca en M Mv Mp

import numpy as np
from scipy import misc
tab = misc.imread('profil.png') # 41*41*4 array

print(tab[20,30,:])   #[y,x,couleur]

#Mv est le profil entier , Mp que l'interieur (sans frontiere)

M=np.zeros([41,41]) #black = 1  white = 0 red=2
Mv=np.zeros([41,41])#black = 1  white = 0 red=1
Mp=np.zeros([41,41])#black = 0  white = 0 red=1

for i in range(41):
    for j in range(41):
        if(tab[i,j,0]==0):
            print('black')
            M[i,j]=1
            Mv[i,j]=1
            Mp[i,j]=0
        if(tab[i,j,0]==255):
            if(tab[i,j,1]==0):
                print('red')
                M[i,j]=2
                Mv[i,j]=1
                Mp[i,j]=1
