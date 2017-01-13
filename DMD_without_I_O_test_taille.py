# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:50:19 2016

@author: l.guillot-salomon
"""

#IMPORTATION DES MODULES

import numpy as np


#CONSTANTES UTILES


timestep = 250*3.5e-9
diameter = 3.78*1.67e-3
speed = 937.7594
nombre_total_fichiers = 256 #Nombre de snapshot d'Arnaud


#LECTURE DE L'INTEGRALITE DU REPERTOIRE / ICI TEST AVEC DES SNAPSHOT RANDOM
  
M = np.random.sample([4998,4]) #On reprend la taille des snapshots d'Arnaud
U1 = np.empty([np.shape(M)[0],nombre_total_fichiers -1])
V1 = np.empty([np.shape(M)[0],nombre_total_fichiers -1])
U2 = np.empty([np.shape(M)[0],nombre_total_fichiers -1])
V2 = np.empty([np.shape(M)[0],nombre_total_fichiers -1])
    
for k in range(nombre_total_fichiers -1): 
    M = np.random.sample([4998,4]) #on crée chaque nouveau snapshot
     
    if k < nombre_total_fichiers -1:
        U1[:,k] = M[:,2]
        V1[:,k] = M[:,3]
            
    if k > 0:
        U2[:,k-1]=M[:,2]
        V2[:,k-1]=M[:,3]
    
    

## calcul sur les champs fluctuants
#U1mean = np.mean(U1)
#V1mean = np.mean(V1)
#U2mean = np.mean(U2)
#V2mean = np.mean(V2)
# 
#for i in range(nombre_total_fichiers -1-1):
#     U1[:,i] = U1[:,i]-U1mean
#     V1[:,i] = U1[:,i]-V1mean
#     U2[:,i] = U2[:,i]-U1mean
#     V2[:,i] = V2[:,i]-U1mean
    

# RESOLUTION DU SYSTEME LINEAIRE
Uop =  np.linalg.lstsq(U1,U2)[0]
Vop =  np.linalg.lstsq(V1,V2)[0]

# CALCUL DES VALEURS ET VECTEURS PROPRES
U_lambda_temp,U_vector = np.linalg.eig(Uop)
V_lambda_temp,V_vector = np.linalg.eig(Vop)

#le programme matlab travaille avec la matrice diagonale des valeurs propres donc nous aussi
U_lambda = np.diag(U_lambda_temp)
V_lambda = np.diag(V_lambda_temp)

# CALCUL DES MODES DMD
Udmd = np.dot(U1,U_vector)
Vdmd = np.dot(V1,V_vector)

# CALCUL DE L'INTENSITE DES MODES NON TRIES (TEMP)
U_mode_energy_temp = np.empty([nombre_total_fichiers -1,1])
V_mode_energy_temp = np.empty([nombre_total_fichiers -1,1])

for i in range(nombre_total_fichiers -1-1):
    U_mode_energy_temp[i] = np.linalg.norm(Udmd[:,i],2)
    V_mode_energy_temp[i] = np.linalg.norm(Vdmd[:,i],2)
 
# CLASSEMENT DES MODES (suivant leur intensit)
 
U_mode_energy = np.sort(U_mode_energy_temp,0)
V_mode_energy = np.sort(V_mode_energy_temp,0)

Uindices = np.argsort(U_mode_energy_temp,0)  
Vindices = np.argsort(V_mode_energy_temp,0) 
 
Utemp = np.empty(np.shape(Udmd),dtype=complex)
Vtemp = np.empty(np.shape(Vdmd),dtype=complex)

for i in range(nombre_total_fichiers -1-1):
    Utemp[:,i] = Udmd[:,Uindices[i,0]]
    Vtemp[:,i] = Vdmd[:,Vindices[i,0]]
    
Udmd = Utemp
Vdmd = Vtemp

Utemp3 = np.empty(np.shape(U_vector),dtype=complex)
Vtemp3 = np.empty(np.shape(V_vector),dtype=complex)

for i in range(nombre_total_fichiers -1-1):
    Utemp3[:,i] = U_vector[:,Uindices[i,0]]
    Vtemp3[:,i] = V_vector[:,Vindices[i,0]]
    
U_vector = Utemp3
V_vector = Vtemp3

# CLASSEMENT DES VALEURS PROPRES
Utemp2 = np.empty(U_lambda.shape, dtype=complex)
Vtemp2 = np.empty(U_lambda.shape, dtype=complex)

for i in range(nombre_total_fichiers -1-1):
    j = Uindices[i,0]
    k = Vindices[i,0]
    Utemp2[i,i] = U_lambda[j,j]
    Vtemp2[i,i] = V_lambda[k,k]
    
U_lambda = Utemp2
V_lambda = Vtemp2

# CALCUL DES FREQUENCES
U_mode_frequency = np.empty([nombre_total_fichiers -1,1])
V_mode_frequency = np.empty([nombre_total_fichiers -1,1])

for i in range(nombre_total_fichiers -1-1):
    
    if np.imag(U_lambda[i,i]) != 0:
        U_mode_frequency[i] = np.imag(np.log2(U_lambda[i,i]))/timestep/2/(np.pi)*diameter/speed
    
    if np.imag(V_lambda[i,i]) != 0:
        V_mode_frequency[i] = np.imag(np.log2(V_lambda[i,i]))/timestep/2/(np.pi)*diameter/speed

 

# CALCUL DES COEFFICIENTS TEMPORELS DES DIFFERENTS MODES
Uphi = np.empty(np.shape(Udmd),dtype=complex)
Vphi = np.empty(np.shape(Vdmd),dtype=complex)

for i in range(nombre_total_fichiers -1-1):
    Uphi[:,i]= Udmd[:,i] / np.linalg.norm(Udmd[:,i],2)
    Vphi[:,i]= Vdmd[:,i] / np.linalg.norm(Vdmd[:,i],2)

Ua=np.dot(np.transpose(U1),Uphi)
Va=np.dot(np.transpose(V1),Vphi)

# SEPARATION PARTIES REELLES ET COMPLEXES
Ureal = np.real(Udmd);
Uimag = np.imag(Udmd);
Vreal = np.real(Vdmd);
Vimag = np.imag(Vdmd);


