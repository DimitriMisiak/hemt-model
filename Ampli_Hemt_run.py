#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 09:35:42 2019

@author: filippini

Gain Hemt
"""

import numpy as np

import matplotlib.pyplot as plt

import time

from scipy.optimize import curve_fit

import matplotlib.patheffects as pe


def fit(t, a, c):
    """
    Sinus function to fit datas
    """
    b = 2 * np.pi / point_periode # b fix the period of sinus
    
    return a * np.sin(b * t + c)


def fit_gauss(x, mu, sigma, ampli):
    """
    Gaussian fit to have an estimation of errors
    """
    return (np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) *
            (sigma * np.sqrt(2 * np.pi))**(-1) * ampli)


plt.close("all")

t = time.time()

# palette couleur Julien###

col = ['black', 'orangered', 'green', 'royalblue', 'fuchsia', 'orange',
       'darkcyan']


# %% INFOS DE BASE SUR LES RUNS


datapath = '/home/filippini/Documents/DATA/RUN59/txt/'

datapath2 = '/home/filippini/Documents/DATA/RUN59'

name_run = 'Evt_20190911_17h01.BIN0_T=4K_fs=10000.0Hz_ts=1s_HEMT_Evt'

name_run2 = 'Evt_20190911_17h01.BIN1_T=4K_fs=10000.0Hz_ts=1s_HEMT_Evt'

info_run = 'test_'

Evtmatrix = np.loadtxt(datapath + name_run)

Evtmatrix2 = np.loadtxt(datapath + name_run2)

# freq sinus
freq = 1

point_periode = 200

size = np.size(Evtmatrix2)

# Debut du fit des donnees avec une periode correspondant a 100 points
periode = 1000
periode2 = 1000
pas = 100

# Fit de la premiere periode visuel p0 a adapter
sinus, pcov = curve_fit(fit, np.arange(0, periode, 1), Evtmatrix[0:periode],
                        sigma = 1.5 + 0 * Evtmatrix[0:periode])

sinus2, pcov2 = curve_fit(fit, np.arange(0, periode, 1), Evtmatrix2[0:periode],
                          sigma = 1.5 + 0 * Evtmatrix[0:periode])
    


def fig():
    """
    Figure du fit avec les datas est la valeur de l'amplification pour
    1 periode
    """
    
    plt.subplot(2,1,1)
    
    # calcul de l'amplification
    
    ampli = sinus2[0] / sinus[0]
    
    # amplification error from curvefit
    
    error_a = pcov[0][0] 
    
    error_a2 = pcov2[0][0] 

    error_amplitude = np.sqrt((error_a2 / sinus[0]**2)  + 
                          (sinus2[0] / sinus[0]**2)**2 * error_a)
    
    print(error_amplitude, np.sqrt((error_a2 / sinus[0]) **2),
          np.sqrt((sinus2[0] * error_a / sinus[0]**2)**2))
    
    # plot signal entree sortie
    plt.plot(Evtmatrix, 'o', color = 'black',
             label='Input signal')
    
    plt.plot(Evtmatrix2, 'o', color='green',
            label='Output signal')
    
    # cut_valeur = np.abs(Evtmatrix) > 5
    
    # amplification = Evtmatrix2[cut_valeur] / Evtmatrix[cut_valeur]
    
    # plt.plot(amplification, linestyle='-', color='gold', linewidth=2,
    #        label='amplification  '+info_run)

    #
    # Plot fit avec data et ampli
    for v in [sinus, sinus2]:
       
        range_sinus = np.arange(0, size, 4)
        
        plt.plot(range_sinus, fit(range_sinus, *v), '-',
                linewidth=2, label='Fit sinus a={0:.3} b={1:.2} a*sin(t+b)'
                .format(*v))
    
    plt.ylabel('Signal [ADU]')
  
    plt.legend()
    
    plt.axis([0,3000, -30, 30])
    
    plt.title('Gain = {0:.3} error = {1:.3}'
             .format(ampli, error_amplitude))
   
    plt.legend(loc='upper right', fontsize='small')


def histo():
    """
    Affiche l'histo de la gaussienne 
    """
    
    # Amplification pour tout les points
    cut_valeur = np.abs(Evtmatrix) > 3

    amplification = (Evtmatrix2[cut_valeur] /
                     Evtmatrix[cut_valeur])

    hist, bin_edges = np.histogram(amplification, 100)
    
    center_bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    fit = curve_fit(fit_gauss, center_bins, hist)
        
    # plot 2
    
    plt.subplot(2,1,2)
    
    # Histogramme de l'amplification
    plt.hist(amplification, 100, label='histogramme')
    
    plt.legend()
    
    range_gaus = np.arange(fit[0][0]-3, fit[0][0]+3, 0.01)
        
    plt.plot(range_gaus,
            fit_gauss(range_gaus, fit[0][0], fit[0][1],
                      fit[0][2]),
            linewidth=2,
            color='darkmagenta',
            path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()],
            label="Fit gaussienne $\\mu$={0:.2e} $\\sigma$={1:2e}".format(*fit[0]))
    
    plt.title('Histogramme')
    
    # plt.axis([-50, 50, 0, 100])
    
    plt.xlabel('Mean')
    
    plt.ylabel('count')
    
    plt.legend(loc='upper right', fontsize='small')


def save():
    
    plt.savefig('/home/filippini/Documents/plot/RUN55/ampli'+info_run+'.png')

plt.figure('amplification', figsize = [20, 10])

fig()

histo()

#save(fig)
