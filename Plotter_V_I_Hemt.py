#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 13:53:47 2019

@author: filippini

Programme qui calcul les VdeI d'un fichier tsv marche pas automatiquement il
faut choisir les lignes

Pas organiser

"""

import numpy as np

import matplotlib.pyplot as pl


datapath = '/home/filippini/Documents/DATA'  # path des datas

runpath = '/home/filippini/Documents'  # path du run

col = ['darkcyan', 'orange', 'green', 'royalblue', 'fuchsia',
       'orange', 'darkcyan']

Vgs = ['100', '170', '150']

col2 = ['green', 'red', 'fuchsia']

symb = ['o', 'v', 's', 'd', '^', 'P', 'X']

# point du spectre de bruit

va_list = [0.141, 0.105, 0.195, 0.111]
vb_list = [0.15, 0.104, 0.199, 0.111]
ib_list = [1.01, 0.98, 1, 0.39]
ia_list = [0.99, 0.98, 1, 0.4]
date_VdeI = ['9 mai 2019', '10 mai 2019']

pl.close("all")

# grille utilisé dans toutes les courbes


def NiceGrid():

    pl.grid(b=True, which='major', color='black', linestyle='-')

    pl.grid(b=True, which='minor', color='silver', linestyle=':')


VdeIpath = datapath + '/RUN55/RUN55_VdeI.tsv'  # path des données

VdeI = np.genfromtxt(VdeIpath, delimiter='\t')


def Fig():
    """
    Figures

    Plot
    ========

    V de I des deux hemts avec des infos concernant la date des prises
    de donnees,
    la valeur des plateaux Isat ainsi que la tension grille source
    et les points des spectres de bruits
    """
    Fig1 = pl.figure("Fig1", figsize=(18 / 1.2, 8 / 1.2))

    Fig1.suptitle('IPNL RUN55 Hemt V(I) voie A et B', fontsize=12,
                  fontweight='bold')

    pl.subplot(121)  # voie A

    j = 0

    axis = [-0, 0.4, 0, 1.2]

    # lignes et colonnes contenant les infos interessantes
    Ids_list = [9, 22]
    Vds_list = [7, 20]

    # donnees voie A pourrait etre fait de facon auto
    Idsmax_list = [1, 0.4]
    Vgsmax_list = [100, 120]

    a = 'nan'  # les casses ne contenant aucune info sont init a nan

    # init
    testids = np.zeros([np.size(Ids_list), int(np.size(VdeI[:][9]) / 2)])
    testvds = np.zeros([np.size(Vds_list), int(np.size(VdeI[:][7]) / 2)])
    j = 0
    maxi = 0
    zero = 0

    # boucle pour creer les tableaux contenant le i de la voie A
    for l in range(np.size(Ids_list)):

        maxi = max(VdeI[Ids_list[l]][:])

        zero = 0
        for i in range(int((np.size(VdeI[:][9])/2))):

            if str(VdeI[Ids_list[l]][j]) != str(a):

                testids[l][i] = VdeI[Ids_list[l]][j]

                if VdeI[Ids_list[l]][j] == 0 and zero > 4:
                    testids[l][i] = maxi
                zero += 1

            j += 2
        j = 0

    j = 0

    maxi = 0
    zero = 0

    # boucle pour creer les tableaux contenant le i de la voie A
    for l in range(np.size(Vds_list)):
        zero = 0
        maxi = max(VdeI[Vds_list[l]][:])

        for i in range(int((np.size(VdeI[:][7]) / 2))):

            if str(VdeI[Vds_list[l]][j]) != str(a):

                testvds[l][i] = VdeI[Vds_list[l]][j]

                if int(VdeI[Vds_list[l]][j]) == 0 and zero > 4:
                    testvds[l][i] = maxi
                zero += 1
            j += 2
        j = 0

    # boucle sur tout les points

    # Plot trop cool voie A
    for i, idsmax, vgsmax, c, date, s in zip(range(np.size(Ids_list)),
                                             Idsmax_list, Vgsmax_list, col,
                                             date_VdeI, symb):

        pl.plot(testvds[i][:]*1e-3, testids[i][:], str(s)+'-', markersize=4,
                color=c, label=date+' Isat='+str(idsmax)+'mA Vgsmax='
                + str(vgsmax)+'mV')

    # Plot des points isole
    pl.plot(va_list, ia_list, "X", markersize=8,
            label='Spectre de bruit 9-10 mai ampli sr560 et ampli sr5184',
            color='black')

    j = 0

    pl.ylabel('Drain source courant [mA]', fontsize=12)
    pl.xlabel('Drain source tension [V]', fontsize=12)
    pl.axis(axis)
    NiceGrid()
    pl.legend()
    pl.legend(loc='lower right', fontsize='x-small')

    pl.subplot(122)
    # plot voie B

    # ini des donnees
    j = 1
    Idsmax_list = [1, 0.4]
    Vgsmax_list = [100, 120]
    NiceGrid()

    testids = np.zeros([np.size(Ids_list), int(np.size(VdeI[:][9]) / 2)])
    testvds = np.zeros([np.size(Vds_list), int(np.size(VdeI[:][7]) / 2)])
    j = 1
    maxi = 0
    zero = 0

    for l in range(np.size(Ids_list)):

        zero = 0
        for i in range(int((np.size(VdeI[:][9]) / 2))):

            if str(VdeI[Ids_list[l]][j]) != str(a):

                maxi = np.max(testids[l][:])
                testids[l][i] = VdeI[Ids_list[l]][j]

                if VdeI[Ids_list[l]][j] == 0 and zero > 4:
                    testids[l][i] = maxi
                zero += 1
            j += 2
        j = 1

    j = 1

    maxi = 0

    for l in range(np.size(Vds_list)):

        zero = 0

        maxi = max(VdeI[Vds_list[l]][:])

        for i in range(int((np.size(VdeI[:][7]) / 2))):

            if str(VdeI[Vds_list[l]][j]) != str(a):

                testvds[l][i] = VdeI[Vds_list[l]][j]

                if VdeI[Vds_list[l]][j] == 0 and zero > 4:
                    testvds[l][i] = maxi
                zero += 1
            j += 2
        j = 1

    for i, idsmax, vgsmax, c, date, s in zip(range(np.size(Ids_list)),
                                             Idsmax_list, Vgsmax_list,
                                             col, date_VdeI, symb):

        pl.plot(testvds[i][:]*1e-3, testids[i][:], str(s)+'-',
                markersize=5, color=c,
                label=date+' Isat='+str(idsmax)+'mA Vgsmax='+str(vgsmax)+'mV')

    pl.plot(vb_list, ib_list, "X", markersize=8,
            label='Spectre de bruit 9-10 mai ampli sr560 et ampli sr5184',
            color='black')

    pl.ylabel('Drain source courant [mA]', fontsize=12)

    pl.xlabel('Drain source tension [V]', fontsize=12)

    pl.axis(axis)

    pl.legend()

    pl.legend(loc='best', fontsize='x-small')

    Fig2 = pl.figure("Fig1", figsize=(18 / 1.2, 8 / 1.2))

    Fig2.suptitle('IPNL RUN54 Hemt V(I) voie A et B 4K',
                  fontsize=12, fontweight='bold')

    return


Fig()
