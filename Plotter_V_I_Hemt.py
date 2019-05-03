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

v_list = [0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1]

v_list2 = [0.3, 0.2, 0.1]

ia_list2 = [3.4, 0.83, 1.46]

ib_list2 = [2.1, 0.69, 1.26]

ib_list = [0.32, 0.95, 1.88, 0.32, 0.93, 1.87, 0.29, 0.78, 1.68]

ia_list = [1.08, 2.3, 3.95, 1.03, 2.3, 3.8, 0.89, 1.93, 3.07]

date_VdeI = ['15 Avril 2019', '15 Avril 2019', '15 Avril 2019',
             '15 Avril 2019', '16 Avril 2019']

pl.close("all")

# grille utilisé dans toutes les courbes


def NiceGrid():

    pl.grid(b=True, which='major', color='black', linestyle='-')

    pl.grid(b=True, which='minor', color='silver', linestyle=':')


VdeIpath = datapath + '/RUN54_VdeI3.tsv'  # path des données

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

    Fig1.suptitle('IPNL RUN54 Hemt V(I) voie A et B', fontsize=12,
                  fontweight='bold')

    pl.subplot(121)  # voie A

    j = 0

    axis = [-0, 2.5, 0, 10]

    # lignes et colonnes contenant les infos interessantes
    Ids_list = [9, 21, 33, 45, 57]
    Vds_list = [7, 19, 31, 43, 55]

    # donnees voie A pourrait etre fait de facon auto
    Idsmax_list = [8.8, 4.3, 4.2, 2.4, 1.2]
    Vgsmax_list = [0, 49, 102, 127, 151]

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
    pl.plot(v_list, ia_list, "X", markersize=8,
            label='Spectre de bruit 17 avril ampli sr560', color='black')

    for v, i, vgs, c in zip(v_list2, ia_list2, Vgs, col2):
        pl.plot(v, i, "P", markersize=8,
                label='Spectre de bruit 19 avril ampli sr5184 Vgs='+str(vgs)
                + 'mV', color=c)

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
    Idsmax_list = [7.2, 4, 2, 1, 0.4]
    Vgsmax_list = [0, 53, 101, 126, 151]
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
        # print(testvds[i][:]*1e-3,testids[i][:])

    pl.plot(v_list, ib_list, "X", markersize=8,
            label='Spectre de bruit 17 avril ampli sr560', color='black')

    for v, i, vgs, c in zip(v_list2, ib_list2, Vgs, col2):

        pl.plot(v, i, "P", markersize=8, label='Spectre de bruit 19' +
                'avril ampli sr5184 Vgs'+str(vgs)+'mV', color=c)

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
