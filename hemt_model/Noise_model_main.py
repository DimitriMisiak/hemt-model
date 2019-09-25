#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 16:34:59 2019

@author: filippini

Programme test qui montre le fonctionnement du code Noise_model_Hemt_ionisa
"""

from Noise_model_Hemt_ionisa import resolution_f
import numpy as np
import matplotlib.pyplot as plt

def NiceGrid():
    """
    Affiche une grille pour un plot,
    evite de marquer à chaque fois ces lignes
    """
    plt.grid(b=True, which='major', color='black', linestyle='-')

    plt.grid(b=True, which='minor', color='silver', linestyle=':')
    
class HEMT:
    """
    class contenant toutes les infos du Hemt et de son bruit
    """
    def __init__(self, e0_ini, ea_ini, eb_ini, i0_ini, ia_ini,
                 ib_ini, Chemt_ini):
        """
        Initialisation des parametres
        """
        self.e0 = e0_ini*1e-9
        self.ea = ea_ini*1e-9
        self.eb = eb_ini*1e-9
        self.i0 = i0_ini*1e-18
        self.ia = ia_ini*1e-18
        self.ib = ib_ini*1e-18
        self.Chemt = Chemt_ini*1e-12

    def in_(self, freq):
        """
        Calcul du bruit en courant du Hemt en fonction de la fréquence
        """
        i_n = np.sqrt(self.i0**2+self.ia**2*freq+(self.ib*freq)**2)
        return i_n

    def en_(self, freq):
        """
        Calcul du bruit en tension du Hemt en fonction de la fréquence
        """
        en = np.sqrt(self.e0**2+self.ea**2/freq+(self.eb/freq)**2)
        return en
# Data list to launch resolution calcul

hemt200 = HEMT(0.18, 5.2, 0, 8.2e-4, 21, 0, 236)
hemt100 = HEMT(0.23, 6.7, 0, 5.3e-4, 16, 0, 103)
hemt40 = HEMT(0.13, 15, 0, 4.7, 7.8, 0, 33)
hemt4 = HEMT(0.22, 36, 0, 4.0e-5, 2.6, 0, 4.6)
hemt2 = HEMT(0.39, 91.4, 0, 3.1, 1.8, 0, 1.8)

# Initialisation des datas

Tb = 20e-3
Rb = 1e10
Cd = 20e-12
Cp = 10e-12
Cc = 2e-9
Cfb = 1e-12

plt.close("all")
# Calcul la résolution en affichant les détails et en rajoutant une source 
# de bruit bias
res = resolution_f(hemt4, Tb, Rb, Cd, Cp, Cc, Cfb, ebias=10e-9, detail = True)

print(res)

# Calcul de la résolution pour différente valeur de la capa det 
for hemt in [hemt2, hemt4, hemt40, hemt100, hemt200]:
    
    res = 0
    
    for Cd in np.arange(1,100,1):
        
        Cd = Cd * 1e-12
        res = np.append(res, resolution_f(hemt, Tb, Rb, Cd, Cp, Cc, Cfb))
    
    res = res[1:]
    
    plt.figure('Scan capa')
    plt.loglog(np.arange(1,100,1), res,
               label = 'Scan hemt {0:.2} pF'.format(hemt.Chemt))
NiceGrid()
plt.xlabel('Cdetector[pF]')
plt.ylabel('Resolution [eV]')
# pl.axis([0, 100, 0, 50])
plt.legend(loc='upper right')
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()