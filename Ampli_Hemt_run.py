#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 09:35:42 2019

@author: filippini

on essaye de fit les sinus pour avoir l'amplification des Hemts
"""

import numpy as np

import matplotlib.pyplot as pl

import time

from scipy.optimize import curve_fit

import matplotlib.patheffects as pe


def fit(t, a, c):
    b = 2 * np.pi / point_periode
    return a * np.sin(b * t + c)


def fit_gauss(x, mu, sigma, ampli):
    return (np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) *
            (sigma * np.sqrt(2 * np.pi))**(-1) * ampli)


pl.close("all")
t = time.time()

# palette couleur Julien###

col = ['black', 'orangered', 'green', 'royalblue', 'fuchsia', 'orange',
       'darkcyan']


# %% INFOS DE BASE SUR LES RUNS


datapath = '/home/filippini/Documents/DATA/RUN55/txt/'

datapath2 = '/home/filippini/Documents/DATA/RUN55'

name_run = 'Evt_20190510_12h53.BIN0_T=4K_fs=10000.0Hz_ts=3s_HEMT_Evt'

name_run2 = 'Evt_20190510_12h53.BIN1_T=4K_fs=10000.0Hz_ts=3s_HEMT_Evt'

info_run = 'Ids=1mA_Vds=200mV_voieB_10kHz'

Evtmatrix = np.loadtxt(datapath + name_run)

Evtmatrix2 = np.loadtxt(datapath + name_run2)

freq = 500
point_periode = 100000/freq

size = np.size(Evtmatrix2)

# Debut du fit des donnees avec une periode correspondant a 100 points
periode = 100
periode2 = 100
pas = 100

# Fit de la premiere periode visuel p0 a adapter
sinus = curve_fit(fit, np.arange(periode), Evtmatrix[0+1000:periode+1000],
                  p0=[17, 2], sigma=0.01 * Evtmatrix[0:periode])

sinus2 = curve_fit(fit, np.arange(periode2), Evtmatrix2[0+1000:periode2+1000],
                   p0=[-190, -40], sigma=0.01*Evtmatrix2[0:periode2])



def fig():
    """
    Figure du fit avec les datas est la valeur de l'amplification pour
    1 periode
    """
    fig = pl.figure("sinus")
    # calcul de l'amplification
    ampli = sinus2[0][0] / sinus[0][0]

    # plot signal entree sortie
    pl.plot(Evtmatrix, linestyle='-', color='black',
            linewidth=2, label='Signal d entree')
    pl.plot(Evtmatrix2, linestyle='-', color='green', linewidth=2,
            label='Signal de sortie '+info_run)
    cut_valeur = np.abs(Evtmatrix) > 5
    amplification = Evtmatrix2[cut_valeur] / Evtmatrix[cut_valeur]
    pl.plot(amplification, linestyle='-', color='gold', linewidth=2,
            label='amplification  '+info_run)

    #
    # Plot fit avec data et ampli
    for v in [sinus, sinus2]:
        range_sinus = np.arange(0, size, 4)
        pl.plot(range_sinus, fit(range_sinus, v[0][0], v[0][1]), 'P',
                markersize=3, label='Fit sinus a={0:.3} b={1:.2} a*sin(t+b)'
                .format(*v[0]))
    # plot 1
    pl.ylabel('Amp [ADU]', fontsize=12)
    pl.axis([0, 200, -350, 350])
    pl.legend()
    pl.title('Mesure amplification par le fit d un sinus gain={0:.3}'
             .format(ampli), fontsize=12)
    pl.legend(loc='upper right', fontsize='small')

    return fig


def histo():

    # Amplification pour tout les points
    cut_valeur = np.abs(Evtmatrix) > 3

    amplification = (Evtmatrix2[cut_valeur] /
                     Evtmatrix[cut_valeur])

    hist, bin_edges = np.histogram(amplification, 100)
    center_bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    fit = curve_fit(fit_gauss, center_bins, hist)
    print(bin_edges, hist)
    # plot 2
    pl.figure("hist")
    # Histogramme de l'amplification
    pl.hist(amplification, 100, label='histogramme')
    pl.legend()
    range_gaus = np.arange(fit[0][0]-10, fit[0][0]+10, 0.01)
    print(fit[0][0])
    pl.plot(range_gaus,
            fit_gauss(range_gaus, fit[0][0], fit[0][1],
                      fit[0][2]),
            linewidth=2,
            color='darkmagenta',
            path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()],
            label="Fit gaussienne $\\mu$={0:.2e} $\\sigma$={1:2e}".format(*fit[0]))
    pl.title('Histogramme', fontsize=12)
    # pl.axis([-50, 50, 0, 100])
    pl.legend(loc='upper right', fontsize='small')


def save(fig):
    pl.savefig('/home/filippini/Documents/plot/RUN55/ampli'+info_run+'.png')


fig = fig()
histo()
# save(fig)
