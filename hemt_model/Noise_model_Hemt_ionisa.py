#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Model Hemt
Created on Tue Apr 23 09:00:09 2019

@author: filippini

Fonction qui permet de calculer les resolutions en voie chaleur et ionisation
"""

import numpy as np

import matplotlib.pyplot as plt


kb = 1.38e-23  # Constante de Boltzmann
e = 1.6e-19  # Charge

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
    
    
def Z_c(freq, C):

    """
    Calcul de l'impédance du condensateur C

    Parameters
    ---------
    freq : float
        Fréquence utilisé pour calculer l'impédance en Hz

    C    : float
        Capacité électrique du condensateur en Farad

    Returns
    ---------
    Z : float
        Impedance du condensateur à la fréquence f complexe
    """
    if C == 'nan':

        a = np.inf*freq

    else:

        a = 1/(2*np.pi*1j*C*freq)
    return a


def ejohnson(T, R):
    """
    Calcul du bruit thermique d'une résistance(Johnson)

    .. math::
       e = \\sqrt{4k_BTR}


    Parameters
    ---------
    freq : float
        Fréquence en Hz

    R    : float
        Résistance en  :math:`\\mathrm{\\Omega}`

    Returns
    ----------
    ej : float
        Bruit johnson de la résistance
    """
    return np.sqrt(4*kb*T*R)

def Z_b(freq, Rb, Cd, Cp, Cc, Cfb, Chemt):
    """
    Calcul de l'impédance vu par le bruit en courant ib
    provenant de la résistance Rb


    Parameters
    ---------
    freq : float
        Fréquence en Hz

    Rb    : float
        Résistance en  :math:`\\mathrm{\\Omega}`

    Cd, Cp, Cc, Chemt : float
        Impédance en F

    Returns
    ----------
    Z : float
        Impédance vu par lel courant ib en A complexe
    """

    Zb = 1/Rb

    Zd = 1/Z_c(freq, Cd)

    Zp = 1/Z_c(freq, Cp)

    Zc = Z_c(freq, Cc)

    Zfb = 1/Z_c(freq, Cfb)

    Zhemt = 1/Z_c(freq, Chemt)

    Z = ((Zc+(Zfb+Zhemt)**(-1))**(-1)+Zb+Zp+Zd)**(-1)

    Z = Z * ((Zfb + Zhemt) ** (-1) / (Zc + (Zfb + Zhemt) ** (-1)))

    return Z


def Z_n(freq, Rb, Cd, Cp, Cc, Cfb, Chemt, detail=None):
    """
    Calcul de l'impédance vu par le bruit en courant in du Hemt


    Parameters
    ---------
    freq : float
        Fréquence en Hz

    Rb    : float
        Résistance en  :math:`\\mathrm{\\Omega}`

    Cd, Cp, Cc, Chemt : float
        Impédance en F

    Returns
    ----------
    Z : float
        Impédance vu par lel courant ib en A complexe
    """
    Zb = 1/Rb

    Zd = 1/Z_c(freq, Cd)

    Zp = 1/Z_c(freq, Cp)

    Zc = Z_c(freq, Cc)

    Zfb = 1/Z_c(freq, Cfb)

    Zhemt = 1/Z_c(freq, Chemt)
    
    if detail is True:
        
        print('\n\n#####Impedance de base#####\n')
        
        for i in [0,9,99]:

            print('frequency =', i+1,
                  ' Zb = {0:e}  Zd = {1:e}'.format(abs(1/Zb), abs(1/Zd[i])),
                  'Zp = {0:e}  \nZc = {1:e}'.format(abs(1/Zp[i]), abs(Zc[i])),
                  'Zfb = {0:e}'.format(abs(1/Zfb[i])),
                  'Zhemt = {0:e}\n'.format(abs(1/Zhemt[i])))
        print('\n')

    Z = ((Zc+(Zb+Zp+Zd)**(-1))**(-1)+Zfb+Zhemt)**(-1)

    return Z


def total_noise(freq, hemt, Tb, Rb, Cd, Cp, Cc, Cfb,
                ebias = None, efb = None,detail = None):
    """
    Parameters
    ==========

    Freq : float

    hemt : Class

    Rb, Tb : float

    Cd,Cp,Cc,Cfb : float

    Return
    ==========
    noise : float
        Bruit du Hemt en linear PSD V/sqrt(Hz)
    """
        
    Zn = Z_n(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt, detail=detail)

    Zfb = Z_n(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)

    Zb = Z_b(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)


    if ebias != None:
        ib = np.sqrt((ejohnson(Tb, Rb)/Rb) ** 2 + (ebias/Rb) ** 2)
        print("test")
        
    else:
        ib =  np.sqrt((ejohnson(Tb, Rb)/Rb) ** 2)
        
        
    if efb != None:
        in_ = np.sqrt(hemt.in_(freq) ** 2 + (efb/np.abs(Zfb)) **2)
        
        
    else:
        in_ = hemt.in_(freq)
    
    
    noise = np.sqrt((in_ * np.abs(Zn)) ** 2
                    + (hemt.en_(freq)) ** 2
                    + (ib * np.abs(Zb)) ** 2)

    if detail is True :
        
        plot_impedance_tot(hemt, freq, Zn, Zb, Zfb, noise, ib)
        
        plot_impedance(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)
        
        plot_noise(freq, hemt, ib, Zb, Zn, noise, ebias , efb)
                
    return noise


def plot_noise(freq, hemt, ib, Zb, Zn, noise, ebias = None, efb = None):
    """
    plot all noise
    """
    
    plt.figure('Noise')
    
    contri_in = hemt.in_(freq) * np.abs(Zn)

    contri_en = hemt.en_(freq)
    
    contri_ib = ib * np.abs(Zb)
    
    plt.loglog(freq, noise, label='bruit total', linewidth=3, color='green')
    
    plt.loglog(freq, np.abs(333 * e * Zb), label='1 keV evt fft',
               color='black', linestyle='--', linewidth=2)
    
    plt.loglog(freq, contri_in, color='green', linestyle=':',
                  label='contribution in')

    plt.loglog(freq, contri_en, color='green', linestyle='--',
              label='contribution en')

    #box_txt(Fignoise, v)
    
    if ebias != None :
        
        plt.loglog(freq, ebias + 0 * freq, color='red', label='Rb noise'+
                   '+ ebias')
        
    else :
        
        plt.loglog(freq, contri_ib, color='purple', label='contribution de Rb')

        
    if efb != None :
        
        plt.loglog(freq, efb + 0 * freq, color='blue', label='feed back noise')
        
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Noise LPSD [V/$\\sqrt{Hz}$]')
    plt.legend(loc='lower left')
    NiceGrid()


def plot_impedance_tot(hemt, freq, Zn, Zb, Zfb, noise, ib):
    
    print('#####Impedance vue par le courant#####\n')

    for i in [0, 9, 99, 999]:
        print('frequency =', i+1, ' Zn = {0:e}'.format(abs(Zn[i])),
              ' Zb = {0:e}  Zfb = {1:e}\n'.format(abs(Zb[i]), abs(Zfb[i])))
        
    print('\n\n')
    print('#####  Bruit  #####\n')
    for i in [1, 10, 100, 1000]:
       
        print('frequency =', i, ' en = {0:e}'.format(hemt.en_(i)),
              ' in = {0:e}'.format(hemt.in_(i) * np.abs(Zn[i-1])),
              ' ib = {0:e}'.format(ib * np.abs(Zb[i-1])),
              ' total_noise = {0:e}'.format(noise[i-1]))
        
    print('\n\n')
        
    plt.figure('Impedance')
    
    plt.suptitle('Impedance')
    
    name = ['Zn tot', 'Zb tot', 'Zfb tot']
    
    for impedance, name in zip([Zn, Zb, Zfb], name):
        
        
        plt.loglog(freq, np.abs(impedance), label = name)
        
        
def plot_impedance(freq, Rb, Cd, Cp, Cc, Cfb, Chemt):
    
    plt.loglog(freq, Rb+0*freq, label='Zb')
    
    NiceGrid()
    
    name = ['Zd', 'Zp', 'Zc', 'Zfb', 'Zhemt']
    
    for impedance, name in zip([Cd, Cp, Cc, Cfb, Chemt], name):
         
        plt.loglog(freq, np.abs(Z_c(freq, impedance)), label = name)
        
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.axis([0.8, 2e4, 1e5, 1e11])
    plt.legend(loc='best')
        
    
def resolution_f(hemt, Tb, Rb, Cd, Cp, Cc, Cfb, ifb=0, ebias=None,
                 i_range=None, df=None, detail=None, efb=None):
    """
    Calcul de la résolution d'un système d'amplification,
    pour un signal discret en frequentielle avec la méthode des trapèzes.

    Parameters
    ---------
    noise_f : np.array
        Bruit fréquentiel (PSD) du module d'amplification en :math:`V^2/Hz`

    i_range : float
        Valeur de la fmax d'integration pour le calcul de la resolution

    Returns
    ----------
    res : float
        resolution du système d'amplification en :math:`eV`
    """
    
    if i_range is not None:
     
        freq = np.arange(i_range[0], i_range[1], df)
    
    if i_range is None:
    
        freq = np.arange(1, 50000, df)
    
    Z = Z_b(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)
    
    noise_f = total_noise(freq, hemt, Tb, Rb, Cd, Cp, Cc, Cfb,
                          ebias=ebias, efb=efb, detail=detail)
    
    noise_f = noise_f**2
    
    Z = np.abs(Z) 
    
    signal_f = 333 * e * Z
    
    if detail is True:
        
        print('#####  FFT  #####\n')
        
        for i in [0,9,99]:

            print('frequency =', i+1,
                  ' FFT = {0:e}'.format(signal_f[i]))
        print('\n')

    fmax = np.size(Z)

    NEPsquare2 = noise_f/(signal_f**2)
    # la case 0 correspond à 1Hz
    NEPsquare2 = 4/NEPsquare2

    # pl.loglog(np.arange(1,200),NEPsquare2,label='2nd methode PSD')

    if i_range is None:

        i_range = [1, fmax]

    if df is None:

        df = 1

    f_min = i_range[0]
    
    f_max = i_range[1]

    res1 = np.sum(NEPsquare2[f_min-1:f_max:df]) ** (-0.5)

    reso = res1*1e3  # On passe en eV
    

    return reso