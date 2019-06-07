#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Model Hemt
Created on Tue Apr 23 09:00:09 2019

@author: filippini

Fonction qui permet de calculer les resolutions en voie chaleur et ionisation
"""

import numpy as np

# import scipy.signal as sgl

from models import model_3exp

import matplotlib.pyplot as pl


kb = 1.38e-23  # Constante de Boltzmann
e = 1.6e-19  # Charge


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


class composant():
    """
    Class contenant les différents para du système
    Test pour voir si c'est pratique
    """

    def __init__(self, Tb_ini, Rb_ini, Cd_ini,
                 Cp_ini, Cc_ini, Cfb_ini):
        """
        Initialisation des para
        """
        self.Tb = Tb_ini
        self.Rb = Rb_ini
        self.Cd = Cd_ini
        self.Cp = Cp_ini
        self.Cc = Cc_ini
        self.Cfb = Cfb_ini
    
    
    def return_all(self):
        """
        Renvoie toutes les valeurs dans un array 
        """
        return [self.Tb, self.Rb, self.Cd, self.Cp, self.Cc, self.Cfb]


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


def Z_fb(freq, Rb, Cd, Cp, Cc, Cfb, Chemt):

    """
    Calcul de l'impédance vu par le bruit en courant ifb


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
        Impédance vu par le courant ib en A
    """

    Zb = 1/Rb

    Zd = 1/Z_c(freq, Cd)

    Zp = 1/Z_c(freq, Cp)

    Zc = Z_c(freq, Cc)

    Zfb = 1/Z_c(freq, Cfb)

    Zhemt = 1/Z_c(freq, Chemt)

    Z = ((Zc+(Zb+Zp+Zd)**(-1))**(-1)+Zfb+Zhemt)**(-1)

    return Z


def Zfet(freq, R, Cfil, Cload):
    """
    Calcul de l'impedance dans le modele de dimitri
    """
    return ((R ** (-1) + Z_c(freq, Cfil) ** (-1) + Z_c(freq, Cload) ** (-1))
            ** (-1))


def total_noise(freq, hemt, Tb=0, Rb=0, Cd=0, Cp=0, Cc=0, Cfb=0,
                ifb=0, detail = None):
    """
    Parameters
    ==========

    Freq : float

    hemt : Class

    Rb, Tb : float

    Cd,Cp,Cc,Cfb : float

    ifb : float

    Return
    ==========
    noise : float
        Bruit du Hemt en linear PSD $V/\\sqrt(Hz)$
    """
    ib = ejohnson(Tb, Rb)/Rb

    Zn = Z_n(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt, detail=detail)

    Zfb = Z_fb(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)

    Zb = Z_b(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)

    noise = np.sqrt((hemt.in_(freq) * np.abs(Zn)) ** 2
                    + (hemt.en_(freq)) ** 2
                    + (ib * np.abs(Zb)) ** 2
                    + (ifb * np.abs(Zfb)) ** 2)



    if detail is True :

        print('#####Impedance vue par le courant#####\n')

        for i in [0, 9, 99]:
            print('frequency =', i+1, ' Zn = {0:e}'.format(abs(Zn[i])),
                  ' Zb = {0:e}  Zfb = {1:e}\n'.format(abs(Zb[i]), abs(Zfb[i])))
            
        print('\n\n')
        print('#####  Bruit  #####\n')
        for i in [1, 10, 100]:
           
            print('frequency =', i, ' en = {0:e}'.format(hemt.en_(i)),
                  ' in = {0:e}'.format(hemt.in_(i) * np.abs(Zn[i-1])),
                  ' ib = {0:e}'.format(ib * np.abs(Zb[i-1])),
                  ' total_noise = {0:e}'.format(noise[i-1]))
            
        print('\n\n')
        
        
    return noise


def pulse_t(fs, detail=None):
    """
    Creation d'un pulse par un modele de 3exp tirer du programme de dimtri
    """
    param = [A1, A2, A3, tau1, tau2, tau3, tauTherm] = [1095.3567632640895,
                                                        444.8384165962467,
                                                        11.646235095672926,
                                                        0.008443239541802526,
                                                        0.03721275993133504,
                                                        663.4475810064974,
                                                        0.008447332918311425]
    
    Gain_adu = 9.8
    param[0] = param[0] * Gain_adu * 1e-9
    param[1] = param[1] * Gain_adu * 1e-9
    param[2] = param[2] * Gain_adu * 1e-9
    t_array = np.arange(0, 1, 1/fs)

    fun = model_3exp(*param, t0=0.5)

    pulse, e1, e2, e3 = fun(t_array, details=True)  # faut plots les trucs

    if detail is True:
        pl.figure("pulse")
        pl.plot(t_array, pulse, 'r')
        pl.plot(t_array, e1, lw=0.5)
        pl.plot(t_array, e2, lw=0.5)
        pl.plot(t_array, e3, lw=0.5)

    return pulse


def resolution_t(noise_f, PSD_signal, i_range=None, fs=None):
    """
    Calcul de la résolution d'un système d'amplification,
    pour un signal discret en temporel.

    Parameters
    ---------
    noise_f : np.array
        Bruit fréquentiel (PSD) du module d'amplification en
         :math:`V^2/Hz`

    signal_t    : np.array
        Signal en temporel d'un pulse de 1keV en :math:`V/keV`

    i_range : float
        Valeur de la frequence max d'integration pour le calcul de
        la resolution

    fs : float
        Frequence d'echantillonnage du signal en Hz

    Return
    ----------
    res : float
        resolution du système d'amplification en :math:`eV`
    """


    pl.figure('res')
    pl.xlabel('Frequency [Hz]', fontsize=22)
    pl.loglog(PSD_signal,label='PSD Run55')
    pl.xticks(fontsize=22)
    pl.yticks(fontsize=22)
    pl.grid(b=True, which='major', color='black', linestyle='-')
    pl.grid(b=True, which='minor', color='silver', linestyle=':')
    pl.axis([1, 200, 1e-24, 1e-18])
    PSD_signal = PSD_signal[1:]
#   fmax = min(np.size(noise_f)-2, int(np.size(signal_t)/2)-1)
    fmax = PSD_signal.shape[0]
    NEPsquare2 = (PSD_signal)/noise_f[1:]

    # pl.loglog(f1, NEPsquare2, label='2nd methode PSD')
    if i_range is None:

        i_range = fmax

    reso2 = 0
    reso2 = np.sum(NEPsquare2[:])
    reso2 = (reso2**(-0.5))
    reso2 = reso2

    return reso2


def resolution_f(noise_f, Z, i_range=None, df=None, detail=None):
    """
    Calcul de la résolution d'un système d'amplification,
    pour un signal discret en frequentielle avec la méthode des trapèzes.

    Parameters
    ---------
    noise_f : np.array
        Bruit fréquentiel (PSD) du module d'amplification en :math:`V^2/Hz`

    Z   : np.array
        Impedance vu par le signal en :math:`\\Omega` complexe

    i_range : float
        Valeur de la fmax d'integration pour le calcul de la resolution

    Returns
    ----------
    res : float
        resolution du système d'amplification en :math:`eV`
    """
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

    res1 = np.sum(NEPsquare2[f_min-1:i_range[1]-1:df]) ** (-0.5)

    reso = res1*1e3  # On passe en eV

    return reso


def res(f_min, f_max, df, hemt, Tb=20e-9, Rb=10e12, Cd=10e-12,
                                               Cp=10e-12, Cc=2e-9, Cfb=1e-12
        , detail= None):
    """
    Calcul de la résolution d'un système d'amplification,
    pour un signal discret en frequentielle.

    Parameters
    ---------
    f_min : np.array
        Borne d'integration inferieur

    f_max   : np.array
        Borne d'integration superieur

    df : float
        Pas d'integration
    Autres para pour modifier les cara du systeme d'amplification
    pour les enelver mettre = 0. Marche pas pour tout.

    Return
    ----------
    res : float
        resolution du système d'amplification en :math:`eV`
    """

    cara = composant(Tb, Rb, Cd, Cp, Cc, Cfb)
    
    freq = np.arange(1, f_max, df)

    compo_list = cara.return_all()
    compo_list = compo_list[1:]

    noise = total_noise(freq, hemt, cara.Tb, *compo_list, detail=detail)

    Z = Z_b(freq, *compo_list, hemt.Chemt)

    reso = resolution_f(noise**2, abs(Z), [f_min, f_max], df, detail=detail)

    return reso
