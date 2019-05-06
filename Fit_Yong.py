#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:46:29 2019

@author: filippini
"""

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


def en_fit(f, e0, ea):
    return np.sqrt(e0 ** 2 + ea ** 2 / f)


def in_fit(freq, ia, i0):
    return np.sqrt(i0 ** 2 + ia ** 2 * freq + (0 * freq) ** 2)


points_inx = np.array([1, 1000])

pointsyong_enx = np.array([1, 10, 100, 1000, 1e9])


# 200pF
points_iny_200 = np.array([21, 680]) * 1e-18


points_eny_200 = np.array([0.52, 0.18]) * 1e-9

points_yong_200 = np.array([5.4, 1.7, 0.52, 0.24, 0.18]) * 1e-9

hemt_200 = '236pF'


# 100pF
points_iny_100 = np.array([15, 510]) * 1e-18

points_alex_iny_100 = np.array([16]) * 1e-18

points_eny_100 = np.array([0.76, 0.22]) * 1e-9

points_yong_100 = np.array([6.3, 2.1, 0.76, 0.34, 0.22]) * 1e-9

para_alex_100 = np.array([0.22, 7.3]) * 1e-9

hemt_100 = '103pF'

# 40pF
points_iny_40 = np.array([9.1, 240]) * 1e-18

points_alex_iny_40 = np.array([9]) * 1e-18

hemt_40 = '36pF'

points_eny_40 = np.array([1.5, 0.12]) * 1e-9

points_yong_40 = np.array([14, 4.5, 1.5, 0.57, 0.12]) * 1e-9

para_alex_40 = np.array([0.12, 16.6]) * 1e-9

# 4pF
points_iny_4 = np.array([2.2, 70]) * 1e-18

points_alex_iny_4 = np.array([2.2]) * 1e-18

points_eny_4 = np.array([4.5, 0.21]) * 1e-9

points_yong_4 = np.array([30, 12, 4.5, 1.4, 0.21]) * 1e-9

para_alex_4 = np.array([0.21, 44.0]) * 1e-9


hemt_4 = '3.6pF'

# 2pF
points_iny_2 = np.array([3.6, 57]) * 1e-18

points_eny_2 = np.array([10, 0.4]) * 1e-9

hemt_2 = '1.8pF'

points_yong_2 = np.array([100, 30, 10, 2.5, 0.4]) * 1e-9

points_enx = np.array([100, 10000000000])

pointseny_list = [points_eny_2, points_eny_4, points_eny_40,
                  points_eny_100, points_eny_200]

pointsiny_list = [points_iny_2, points_iny_4, points_iny_40,
                  points_iny_100, points_iny_200]

pointsyong_list = [points_yong_2, points_yong_4, points_yong_40,
                   points_yong_100, points_yong_200]

paraalex_list = [para_alex_4, para_alex_40, para_alex_100]

paraalexin_list = [points_alex_iny_4, points_alex_iny_40, points_alex_iny_100]

capa_list = [hemt_2, hemt_4, hemt_40, hemt_100, hemt_200]

col = ['black', 'orangered', 'green', 'royalblue', 'fuchsia',
       'orange', 'darkcyan']

test_en_list = [curve_fit(en_fit, pointsyong_enx, v, sigma=0.1 * v)[0]
                for v in pointsyong_list]  # fit de en avec les donnees de Yong

# test_en_list = [curve_fit(en_fit, pointsyong_enx, v, sigma = 0.01)[0]
#                 for v in pointsyong_list]

test_in_list = [curve_fit(in_fit, points_inx, v, sigma=0.1 * v)[0]
                for v in pointsiny_list]  # fit de in avec les donnees de Yong

t_range_full = np.linspace(1, 50000, 50000)

axiss = [0.6, 1e5, 1e-10, 200e-9]

plt.close("all")

figen = plt.figure("en")

plt.title("Fit_en", fontsize=12)

# figen.suptitle('Noise en HEMT', fontsize=12, fontweight = 'bold')

for v, capa, yong, colo in zip(test_en_list, capa_list, pointsyong_list, col):

    v = abs(v)

    plt.loglog(t_range_full, en_fit(t_range_full, *v), color=colo,
               linestyle='-', label='Hemt C=' + capa +
                                    ' e0={0:.2e}, ea={1:.2e}'.format(*v))

    plt.plot(pointsyong_enx, yong, linestyle='', color=colo, marker='o',
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

plt.close()
figin = plt.figure("in")

# figin.suptitle('Noise in HEMT', fontsize=12, fontweight = 'bold')
plt.title("Fit_in", fontsize=12)

for v, capa, iny, colo in zip(test_in_list, capa_list, pointsiny_list, col):

    v = abs(v)
    plt.loglog(t_range_full, in_fit(t_range_full, *v),  linestyle='-',
               color=colo, label='Hemt C=' + capa
                                 + '  i0={1:.2e}, ia={0:.2e}'.format(*v))

    plt.plot(points_inx, iny, color=colo, linestyle='', marker='o',
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
