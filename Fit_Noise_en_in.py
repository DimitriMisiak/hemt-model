#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:46:29 2019

@author: filippini

Fit des datas pour modeliser le bruit des hemts
"""

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


def en_fit(f, e0, ea):
    """
    Return le bruit en tension "en" fonction des para e0 ea et la freq
    """
    return np.sqrt(e0 ** 2 + ea ** 2 / f)


def in_fit(freq, ia, i0):
    """
    Return le bruit en courant "in" fonction des para i0 ia et la freq
    """
    return np.sqrt(i0 ** 2 + ia ** 2 * freq + (0 * freq) ** 2)


# Points yong article pour le fit

points_en_yong_article_92 = np.array([16, 5.3, 1.8, 0.6]) * 1e-9
points_in_yong_article_92 = np.array([15, 49, 160, 510]) * 1e-18

points_en_yong_article_26 = np.array([26, 8, 2.5, 0.8]) * 1e-9
points_in_yong_article_26 = np.array([9.1, 27, 80, 240]) * 1e-18

points_en_yong_article_5 = np.array([40, 14, 4.8, 1.7]) * 1e-9
points_in_yong_article_5 = np.array([2.8, 9.2, 28, 100]) * 1e-18

#freq pour le fit
point_yong_article = np.array([1, 10, 100, 1000])

point_list_article_en = [points_en_yong_article_92, points_en_yong_article_26,
                         points_en_yong_article_5]

point_list_article_in = [points_in_yong_article_92, points_in_yong_article_26,
                         points_in_yong_article_5]

capa_list = ['92pF', '26pF', '5.3pF']

col = ['black', 'orangered', 'green', 'royalblue', 'fuchsia',
       'orange', 'darkcyan']


# Fit du bruit 
test_en_list = [curve_fit(en_fit, point_yong_article, v, sigma=0.01 * v)[0]
                for v in point_list_article_en]

test_in_list = [curve_fit(in_fit, point_yong_article, v, sigma=0.01 * v)[0]
                for v in point_list_article_in]



t_range_full = np.linspace(1, 50000, 50000)
axiss = [0.6, 1e5, 1e-10, 200e-9]
plt.close("all")
figen = plt.figure("en")
plt.title("Fit_en", fontsize=12)
# figen.suptitle('Noise en HEMT', fontsize=12, fontweight = 'bold')

# Plot le fit en et les para
for v, capa, yong, colo in zip(test_en_list, capa_list, point_list_article_en,
                               col):

    v = abs(v)

    plt.loglog(t_range_full, en_fit(t_range_full, *v), color=colo,
               linestyle='-', label='Hemt C=' + capa +
                                    ' e0={0:.2e}, ea={1:.2e}'.format(*v))

    plt.plot(point_yong_article, yong, color=colo, marker='o', linestyle=' ',
             markersize=5, label='Donnees en table 1 Hemt de '+capa)

#    plt.loglog(t_range_full, en_fit(t_range_full, *para), color=colo,
#               linestyle=':', label='Fit_en Alex Hemt C=' + capa +
#                                  ' e0={0:.3e} V, ea={1:.3e} V'.format(*para)
    plt.xlabel('Freq [Hz]', fontsize=12)

    plt.ylabel('Noise [V/$\\sqrt{Hz}$]', fontsize=12)


plt.grid(b=True, which='major', color='black', linestyle='-')

plt.grid(b=True, which='minor', color='silver', linestyle=':')

plt.legend()

plt.legend(loc='upper right', fontsize='x-small')
plt.axis(axiss)


# print('e0=',test_en[0][0],' ea=',test_en[0][1])

figin = plt.figure("in")

# figin.suptitle('Noise in HEMT', fontsize=12, fontweight = 'bold')
plt.title("Fit_in", fontsize=12)


# Plot le fit en et les para
for v, capa, iny, colo in zip(test_in_list, capa_list, point_list_article_in,
                              col):

    v = abs(v)
    plt.loglog(t_range_full, in_fit(t_range_full, *v),  linestyle='-',
               color=colo, label='Hemt C=' + capa
                                 + '  i0={1:.2e}, ia={0:.2e}'.format(*v))

    plt.plot(point_yong_article, iny, color=colo, linestyle='', marker='o',
             markersize=5, label='Donnees in table 1 Hemt de '+capa)

    # plt.loglog(t_range_full, in_fit(t_range_full,*para),color=colo,
    #           linestyle=':', label='Fit_in  Alex Hemt C=' + capa +
    #                           '     ia={0:.3e}'.format(*para))

    plt.xlabel('Freq [Hz]', fontsize=12)

    plt.ylabel('Noise [A/$\\sqrt{Hz}$]', fontsize=12)


plt.grid(b=True, which='major', color='black', linestyle='-')
plt.grid(b=True, which='minor', color='silver', linestyle=':')

plt.legend()

plt.legend(loc='upper left', fontsize='x-small')
