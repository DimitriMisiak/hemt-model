#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:08:18 2019

@author: filippini
"""

import numpy as np

import matplotlib.pyplot as pl


class HEMT:

    def __init__(self, e0_ini, ea_ini, i0_ini, ia_ini, ib_ini, Chemt_ini):
        self.e0 = e0_ini * 1e-9
        self.ea = ea_ini * 1e-9
        self.i0 = i0_ini * 1e-18
        self.ia = ia_ini * 1e-18
        self.ib = ib_ini * 1e-18
        self.Chemt = Chemt_ini * 1e-12

    def in1(self, freq):
        self.i_n = np.sqrt(self.i0 ** 2 + self.ia ** 2 * freq
                           + (self.ib * freq) ** 2)
        return self.i_n

    def en1(self, freq):
        self.en = np.sqrt(self.e0 ** 2 + self.ea ** 2 / freq)
        return self.en


class FET:

    def __init__(self, e0_ini, ea_ini, eb_ini, i0_ini, ia_ini,
                 ib_ini, Cfet_ini):
        self.e0 = e0_ini * 1e-9
        self.ea = ea_ini * 1e-9
        self.eb = eb_ini * 1e-9
        self.i0 = i0_ini * 1e-18
        self.ia = ia_ini * 1e-18
        self.ib = ib_ini * 1e-18
        self.Cfet = Cfet_ini * 1e-12

    def in1(self, freq):
        self.i_n = np.sqrt(self.i0 ** 2 + self.ia ** 2 * freq
                           + (self.ib * freq) ** 2)
        return self.i_n

    def en1(self, freq):
        self.en = np.sqrt(self.e0 ** 2 + self.ea ** 2 / freq
                          + (self.eb / freq) ** 2)
        return self.en


# impedance du syteme avant le HEMT

def Zc(freq, C):

    return 1 / (2 * 1j * np.pi * freq * C)


def Z(freq, R, C):

    return abs(1 / ((1 / R) + (1 / Zc(freq, C))))


def ejohnson(freq, T, R, C):

    return abs((Zc(freq, C) / (Zc(freq, C) + R)) * (4 * kb * T * R))


t200 = HEMT(0.18, 5.2, 4.5e-5, 21, 0, 236)

t100 = HEMT(0.23, 6.2, 1.7e-4, 15.5, 0, 103)

t40 = HEMT(0.12, 14.9, 5.0, 7.59, 0, 33)

t4 = HEMT(0.21, 37.1, 4.3e-6, 2.21, 0, 4.6)

t2 = HEMT(0.4, 91.4, 3.1, 1.8, 0, 1.8)

Fet_dimitri = FET(1.6, 23.9, 7.81, 6.17, 0.276, 1.28, 48)

Hemt_list = [t2, t4, t40, t100, t200]

kb = 1.38e-23

col = ['black', 'orangered', 'royalblue', 'fuchsia', 'orange', 'darkcyan']

fs = range(50000)

R = 88e6
T = 13.5e-3


def NiceGrid():

    pl.grid(b=True, which='major', color='black', linestyle='-')

    pl.grid(b=True, which='minor', color='silver', linestyle=':')


def fig1():

    Fig = pl.figure('Fig', figsize=(12.8 / 1.2, 8 / 1.2))

    Fig.suptitle('Hemt ', fontsize=12, fontweight='bold')

    NiceGrid()

    pl.ylabel('Noise [V/$\\sqrt{Hz}$]', fontsize=12)

    # initialisation des données
    en = np.zeros(int(np.size(fs)))

    entot = np.zeros(int(np.size(fs)))

    in_ = np.zeros(int(np.size(fs)))

    impedance = np.zeros(int(np.size(fs)))

    bjohnson = np.zeros(int(np.size(fs)))

    tes = 0
    # calcul du bruit de en et in
    for v in [t100]:
        # while i2<np.size(np.array(R[0][:])):
        for l in range(1, int(np.size(fs) + 1)):

            impedance[l-1] = Z(l, R, v.Chemt)

            en[l-1] = v.en1(l)

            in_[l-1] = v.in1(l)

            bjohnson[l-1] = ejohnson(l, T, R, v.Chemt)

            entot[l-1] = np.sqrt(en[l-1] ** 2 +
                                 (in_[l-1]*impedance[l-1])**2+bjohnson[l-1])

            in_[l-1] *= impedance[l-1]

            '''
            if l in [1,1000]:
                zae=1/abs((1/(R[0][i2]*1e6))+(1/(1j*v.Chemt*l*2*np.pi)))
                print(zae)
                print(str(l)+'Hz Z='+str(impedance[l-1]*1e-6)+'MO   '
                + str(v.Chemt*1e12),' in=',str(a))
                 '''

        pl.loglog(range(1, int(np.size(fs))), entot[1:], linestyle='-',
                  linewidth=2,
                  color=col[tes], label='Hemt pour R=' + str(R * 1e-6) +
                                        'MO C=' + str(v.Chemt * 1e12) + 'pF')

        pl.loglog(range(1, int(np.size(fs))), in_[1:], linestyle=':',
                  color=col[tes], label='in*Z Hemt R=' + str(R * 1e-6) +
                                        'MO C=' + str(v.Chemt*1e12) + 'pF')

        pl.loglog(range(1, int(np.size(fs))), en[1:], linestyle='--',
                  color=col[tes], label='en Hemt R=' + str(R * 1e-6) +
                                        'MO C=' + str(v.Chemt * 1e12)+'pF')

        # pl.loglog(range(int(np.size(fs))), bjohnson, linestyle='-.'
        #          ,color=col[tes],label='')

        tes += 1
    v = Fet_dimitri
    for l in range(1, int(np.size(fs) + 1)):

        impedance[l-1] = Z(l, R*1e6, v.Cfet)

        en[l-1] = v.en1(l)

        in_[l-1] = v.in1(l)

        entot[l-1] = np.sqrt(en[l-1] ** 2 + (in_[l-1] * impedance[l-1]) ** 2)

        in_[l-1] *= impedance[l-1]
        # i2+=2
    pl.loglog(range(1, int(np.size(fs))), entot[1:], linewidth=3,
              linestyle=':', color=col[tes], label='Fet pour R=' + str(R*1e-6)
              + 'MO C=' + str(v.Cfet*1e12) + 'pF')

    pl.legend()
    pl.legend(loc='upper right', fontsize='x-small')


fig1()
