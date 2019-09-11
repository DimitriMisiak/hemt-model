#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 16:33:12 2019

@author: filippini
"""

import numpy as np

import matplotlib.pyplot as plt

from models import model_3exp

import scipy.signal as sgl



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
    

def Zc(freq, C):
    """
    Electrical impedance of capacitor
    """
    return 1/(2*1j*np.pi*freq*C)


def ejohnson(T, R):
    """
    Johnson Noise of resistance
    """
    return (4*kb*T*R)


def total_noise(freq, hemt, T, Rntd, edaq = 0, detail = None):
    """
    Return noise Linear PSD 
    """
    
    Zntd = Rntd
        
    Zhemt = Zc(freq, hemt.Chemt)
    
    entd = np.sqrt(ejohnson(T, Rntd))
    
    intd =  entd / Rntd
    
    Z = Zntd * Zhemt / (Zntd + Zhemt)
    
    noise = np.sqrt( (intd ** 2 + hemt.in_(freq) ** 2) * np.abs(Z) ** 2
                    + (hemt.en_(freq)) ** 2
                    + (edaq) ** 2)
    
    if detail is True:
        
        plot_noise(freq, noise, hemt, edaq, Z, entd *np.abs( Zhemt / (Zntd + Zhemt)))
    
    return noise


def plot_noise(freq, noise, hemt, edaq, Z, entd):
    """
    Noise Plot
    """
    
    plt.figure('Noise')
    
    contri_in = hemt.in_(freq)*np.abs(Z)

    contri_en = hemt.en_(freq)
    
    plt.loglog(freq, noise, label='bruit total', linewidth=3, color='green')
    
    plt.loglog(freq, contri_in, color='green', linestyle=':',
                  label='contribution in')

    plt.loglog(freq, contri_en, color='green', linestyle='--',
              label='contribution en')
    
    plt.loglog(freq, entd + 0*freq, color='blue', label='contribution NTD')
    #box_txt(Fignoise, v)

    plt.loglog(freq, (edaq + 0*freq), color='purple', label='contribution Daq',
               linestyle = '-.')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Noise LPSD [V/$\\sqrt{Hz}$]')
    plt.legend(loc='lower left')
    NiceGrid()
    
    
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
                                                        0.008447332918311425,
                                                        ]
    


    Gain_adu = 9.8 / 5.89
    param[0] = param[0] / Gain_adu * 1e-9
    param[1] = param[1] / Gain_adu * 1e-9
    param[2] = param[2] / Gain_adu * 1e-9
    t_array = np.arange(0, 1, 1/fs)

    fun = model_3exp(*param, t0=0.5)

    pulse, e1, e2, e3 = fun(t_array, details=True)  # faut plots les trucs
    if detail is True:
        plt.figure("pulse")
        plt.plot(t_array, pulse, 'r')
        plt.plot(t_array, e1, lw=0.5)
        plt.plot(t_array, e2, lw=0.5)
        plt.plot(t_array, e3, lw=0.5)
        plt.xlabel('time', fontsize=14)
        plt.ylabel('Amplitude', fontsize=14)
    
    return pulse
        
        
def resolution_calcul(freq, noise_f, signal_t, detail = None):
    
    fs = np.max(freq)
    
    fs = 400
    
    f, signal_f = sgl.welch(signal_t, fs, 'boxcar', int(fs), noverlap=0)
    
    signal_f = signal_f[1:] 
        
    NEPsquare = (signal_f[0:int(fs/2)-1] / (noise_f[0:int(fs/2)-1])**2) 
    
    resolution = np.sum(NEPsquare) ** (-0.5)
    
    #resolution = 1
    
    if detail is True:
        
        #print(signal_f, np.size(signal_f), f)
        plt.figure('signal welch')
        plt.loglog(f[1:int(fs/2)],signal_f[:int(fs/2)-1])
        
    return resolution

def test():
    
    freq = np.arange(1, 50000)
    
    R = 2.25e6
    
    T = 15e-3
    
    hemt200 = HEMT(0.18, 5.2, 0, 4.5e-5, 21, 0, 236)
    
    noise = total_noise(freq, hemt200, T, R, detail = True)
    
    pulse = pulse_t(50000)
    
    res = resolution_calcul(freq, noise, pulse, detail = True)
    
    print('la resolution est de', res, 'eV')
    
    