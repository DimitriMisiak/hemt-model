#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 09:00:09 2019

@author: filippini
"""

import numpy as np

import scipy.signal as sgl


kb = 1.38e-23  # Constante de Boltzmann
e = 1.6e-19  # Charge de l'électron


class HEMT:
    """
    class contenant toutes les infos du Hemt et de son bruit
    """
    def __init__(self, e0_ini, ea_ini, eb_ini, i0_ini, ia_ini,
                 ib_ini, Chemt_ini):

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

    def __init__(self, Tb_ini=0.02, Rb_ini=1e10, Cd_ini=1e-11,
                 Cp_ini=1e-12, Cc_ini=2e-9, Cfb_ini=1e-12):

        self.Tb = Tb_ini
        self.Rb = Rb_ini
        self.Cd = Cd_ini
        self.Cp = Cp_ini
        self.Cc = Cc_ini
        self.Cfb = Cfb_ini


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

        return np.inf*freq

    else:

        return 1/(2*np.pi*1j*C*freq)


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


def Z_b(freq, Rb, Cd=np.inf, Cp=np.inf, Cc=1, Cfb=1, Chemt=1):
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


def Z_n(freq, Rb, Cd, Cp, Cc, Cfb, Chemt):
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


def total_noise(freq, hemt, Tb=0, Rb=0, Cd=0, Cp=0, Cc=0, Cfb=0,
                ifb=0):
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

    Zc = Z_c(freq, Cc)

    Zfb1 = 1/Z_c(freq, Cfb)

    Zhemt = 1/Z_c(freq, hemt.Chemt)

    Zn = Z_n(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)

    Zfb = Z_fb(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)

    Zb = Z_b(freq, Rb, Cd, Cp, Cc, Cfb, hemt.Chemt)

    # print('Zb = ', abs(Zb))
    # print('ib= ', np.abs(Zb)*ib, ' \nen= ', hemt.en_(freq),
    #      ' \nin= ', hemt.in_(freq)*np.abs(Zn))

    noise = np.sqrt((hemt.in_(freq) * np.abs(Zn)) ** 2
                    + (hemt.en_(freq)) ** 2
                    + (ib * np.abs(Zb)) ** 2
                    + (ifb * np.abs(Zfb)) ** 2)

    return noise


def resolution(noise_f, signal_t, i_range=None, fs=None):
    """
    Calcul de la résolution d'un système d'amplification,
    pour un signal discret en temporel avec la méthode des trapèze.

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
    if fs is None:

        fs = np.size(signal_t)

    f1, PSD_signal = sgl.welch(signal_t, fs, 'boxcar', int(1*fs))

    fmax = min(np.size(noise_f)-2, int(np.size(signal_t)/2)-1)

    NEPsquare2 = noise_f[1:fmax]/(PSD_signal[1:fmax])

    NEPsquare2 = 4/NEPsquare2

    # pl.loglog(np.arange(1,200),NEPsquare2,label='2nd methode PSD')

    if i_range is None:

        i_range = fmax

    res = 1/(np.trapz(NEPsquare2, x=np.arange(1, i_range))**(0.5))

    res = res*1e3  # On passe en eV

    return res


def resolution_f(noise_f, Z, i_range=None, df=None):
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
    signal_f = 333*e*Z

    fmax = np.size(Z)

    NEPsquare2 = noise_f/(signal_f**2)
    # la case 0 correspond à 1Hz
    NEPsquare2 = 4/NEPsquare2

    # pl.loglog(np.arange(1,200),NEPsquare2,label='2nd methode PSD')

    if i_range is None:

        i_range = fmax

    if df is None:

        df = 1

    f_min = df

    # print(np.size(NEPsquare2[f_min:i_range]),np.size(freq_trapz))

    res1 = 1/(np.sum(NEPsquare2[f_min-1:i_range-1:df])**(0.5))

    res = res1*1e3  # On passe en eV

    return res


def res(f_min, f_max, df, hemt, composant=composant()):
    """
    Calcul de la résolution d'un système d'amplification,
    pour un signal discret en frequentielle avec la méthode des trapèzes.

    Parameters
    ---------
    f_min : np.array
        Bruit fréquentiel (PSD) du module d'amplification en :math:`V^2/Hz`

    f_max   : np.array
        Impedance vu par le signal en :math:`\\Omega` complexe

    df : float
        Valeur de la fmax d'integration pour le calcul de la resolution

    Returns
    ----------
    res : float
        resolution du système d'amplification en :math:`eV`
    """
    freq = np.arange(1, 50000, 1)
    Rb = composant.Rb
    Tb = composant.Tb
    Cc = composant.Cc
    Cd = composant.Cd
    Cp = composant.Cp
    Cfb = composant.Cfb
    compo_list = [Rb, Cd, Cp, Cc, Cfb]
    noise = total_noise(freq, hemt, Tb, *compo_list, 0)

    Z = Z_b(freq, *compo_list, hemt.Chemt)
    print(abs(Z))
    res = resolution_f(noise**2, abs(Z), f_max, df)
    print(res)
    return res
