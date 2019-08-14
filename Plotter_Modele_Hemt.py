#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 08:54:09 2019

Plot Modele des hemts

@author: filippini
"""

import numpy as np
import matplotlib.pyplot as pl
# import sys
# from os import path
import Modele_ampli_Hemt_bruit as mbh

# sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
# initialisation des trucs pratiques

col = pl.cm.tab10(np.linspace(0, 1, 10))
axis = [1, 1e4, 1e-10, 1e-7]

pl.close("all")

# toutes les def


def NiceGrid():
    """
    Affiche une grille pour un plot,
    evite de marquer à chaque fois ces lignes
    """
    pl.grid(b=True, which='major', color='black', linestyle='-')

    pl.grid(b=True, which='minor', color='silver', linestyle=':')


def Gridplot():
    pl.grid(b=True, which='major', color='black', lw=0.1, linestyle=':')


def box_txt(ax, hemt, loc=[0.70, 0.85]):
    """
    Texte decrivant les donnees du Hemt
    Ainsi que les impedances
    """
    text = ('IPNL-C2N {:.2}F Hemt model with:\n'.format(hemt.Chemt)
            + '$e_n= \\sqrt{e_0^2+e_a^2/f+e_b^2/f^2}$ [$V /\\sqrt{Hz}$]\n'
            + '$i_n=\\sqrt{i_0^2+i_a^2 \\cdot f+i_b^2 \\cdot f^2}$'
            + '[$A/\\sqrt{Hz}$] \n'
            + '$i_0$={:.2}  , $i_a$={:.2} , $i_b$={:.2}\n'
            .format(hemt.i0, hemt.ia, hemt.ib)
            + '$e_0$={:.2}  , $e_a$={:.2} , $e_b$={:.2}\n'
            .format(hemt.e0, hemt.ea, hemt.eb))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

    # place a text box in upper left in axes coords
    ax.text(loc[0], loc[1], text, fontsize=10,
            verticalalignment='top', bbox=props)
    ax.suptitle('Tb={0:}K Rb={1:.2e}$\\Omega$  Cd={2:}F'
                .format(*compo)
                + ' Cp={3:}F Cc={4:}F Cfb={5:}F'
                .format(*compo),
                y=0.92, fontsize=10)


def figure(hemt):
    """
    Figure montrant toutes les contributions pour le bruit du Hemt

    """
    for v in hemt:
        i = 0

        # Calcul du bruit total

        noise = mbh.total_noise(freq, v, Tb, Rb, Cd, Cp, Cc, Cfb, ifb)

        # Creation de la figure
        Fignoise = pl.subplots()
        pl.title('Noise', y=1.07, fontsize=16, fontweight='bold')

        # Plot du bruit total
        pl.loglog(freq, noise, label='bruit total', linewidth=3, color='green')
        i += 1

        # Calcul pour les dif contri
        Zb = (mbh.Z_b(freq, Rb, Cd, Cp, Cc, Cfb, v.Chemt))

        contri_in = v.in_(freq)*np.abs(mbh.Z_n(freq, Rb,
                                               Cd, Cp, Cc, Cfb, v.Chemt))

        contri_en = v.en_(freq)

        res_f = mbh.resolution_f(noise**2, np.abs(Zb))

        # Plot du signal en frequentielle
        pl.loglog(freq, 333*1.6e-19*np.abs(Zb), color='black', linestyle='--',
                  linewidth=2, label='test 1keV')

        contri_ib = np.abs(Zb) * mbh.ejohnson(Tb, Rb) / Rb

        fft = 333*1.6e-19*np.abs(Zb)

        print(res_f)

        pl.loglog(freq, contri_in, color='green', linestyle=':',
                  label='contribution in')

        pl.loglog(freq, contri_en, color='green', linestyle='--',
                  label='contribution en \n res={:.3}eV '.format(res_f))

        #box_txt(Fignoise, v)

        pl.loglog(freq, contri_ib, color='purple', label='contribution de Rb')
        pl.xlabel('Frequency [Hz]')
        pl.ylabel('Noise LPSD [V/$\\sqrt{Hz}$]')
        NiceGrid()

        for i in [1, 10, 100, 1000]:

            print('freq= {0} Hz Chemt= {1}F in= {2:.5} V/Hz en = {3:.5}'
                  .format(i, v.Chemt*1e12, np.abs(mbh.Z_n(i, Rb,
                                               Cd, Cp, Cc, Cfb, v.Chemt)), v.en_(i))
                  + ' fft = {0:.5}'.format(fft[i-1]))

        pl.axis([0.9, 2e4, 1e-10, 1e-7])
        pl.legend(loc='lower left')
        figManager = pl.get_current_fig_manager()
        figManager.window.showMaximized()
        pl.show()


def figure_impedance(hemt):
    """
    Plot les impédances des capa et des resistances

    Parameter
    =========
    hemt : list
        tracer des impedances pour tout les hemts
    """
    for v in hemt:

        Zb_tot = np.abs(mbh.Z_b(freq, Rb, Cd, Cp, Cc, Cfb, v.Chemt))

        Zfb_tot = np.abs(mbh.Z_fb(freq, Rb, Cd, Cp, Cc, Cfb, v.Chemt))

        Zin_tot = np.abs(mbh.Z_n(freq, Rb, Cd, Cp, Cc, Cfb, v.Chemt))

        Fignoise = pl.figure('Impedance '+str(v.Chemt))

        pl.title('Impedance', y=1.07, fontsize=16, fontweight='bold')

        Fignoise.suptitle('Impedance', fontsize=16, fontweight='bold')

        pl.loglog(freq, Zb_tot, label='Zb tot', color='green', linewidth=1)

        pl.loglog(freq, Zfb_tot, label='Zfb tot', color='orange', linewidth=1)

        pl.loglog(freq, Zin_tot, label='Zn tot', color='blue', linewidth=1)

        pl.loglog(freq, Rb+0*freq, label='Zb', linewidth=2)

        pl.loglog(freq, abs(mbh.Z_c(freq, Cd)), label='Zd', linewidth=2)

        pl.loglog(freq, abs(mbh.Z_c(freq, Cp)), label='Zp', linewidth=2)

        pl.loglog(freq, abs(mbh.Z_c(freq, Cc)), label='Zc', linewidth=2)

        pl.loglog(freq, abs(mbh.Z_c(freq, Cfb)), label='Zfb', linewidth=2)

        pl.loglog(freq, abs(mbh.Z_c(freq, v.Chemt)), label='Zhemt',
                  linewidth=2)

        # Para du plot
        pl.xlabel('Frequency [Hz]')
        pl.ylabel('Impedance [$\\Omega$]')
        pl.axis([0.8, 2e4, 1e5, 1e11])
        box_txt(Fignoise, v)
        NiceGrid()
        pl.legend(loc='lower left')
        figManager = pl.get_current_fig_manager()
        figManager.window.showMaximized()
        pl.show()


def fig_scan_Cdet(hemt):
    """
    Figure montrant l'évolution de la résolution en fonction de la capa
    du detecteur

    """
    for v in hemt:

        i = 0

        Fignoise = pl.figure('Scan '+str(v.Chemt))

        pl.title('Impedance', y=1.07, fontsize=16, fontweight='bold')
        for df in [1, 2]:

            max_Cdet = 300

            res = np.zeros(max_Cdet)

            for Cd in np.arange(1, max_Cdet+1, 1):

                j = Cd

                res[j-1] = mbh.res(df, 50000, df, v, Cd=Cd*1e-12)

            pl.loglog(np.arange(1, max_Cdet+1, 1), res, label='Res eV with df=' +
                    str(df)+'Hz', linewidth=2, color=col[i])

            i += 1

            Cd = 'Scan'
            compo[2] = Cd
            box_txt(Fignoise, v, loc=[0.17, 0.85])

            pl.xlabel('Cdetector[pF]')
            pl.ylabel('Resolution [eV]')
            NiceGrid()
            pl.legend(loc='lower right')
            figManager = pl.get_current_fig_manager()
            figManager.window.showMaximized()
            pl.show()


def fig_scan_Cdet_compa(hemt):
    """
    Figure montrant l'évolution de la résolution en fonction de la capa
    du detecteur

    """
    Cd = 'Scan'
    # composant[2] = Cd
    pl.figure('Scan compa')
    '''
    pl.suptitle('Tb={0:}K Rb={1:.2e}$\\Omega$  Cd={2:}'
                .format(*composant)
                + ' Cp={3:}F Cc={4:}F Cfb={5:}F'
                .format(*composant),
                y=0.92, fontsize=10)
    pl.title('Scan Cd', y=1.07, fontsize=16, fontweight='bold')
    '''
    i = 1
    k = 0
    composan = compo.return_all()
    
    print(Cp, Cc)
    for v in hemt:

        if v != Fet_dimitri:
            df = 1

            max_Cdet = 150

            res = np.zeros(max_Cdet)

            for Cd in np.arange(1, max_Cdet+1, 1):

                j = Cd
                
                composan[2] = j * 1e-12
                res[j-1] = mbh.res(df, 50000, df, v, *composan)
                if Cd == 10:
                    print(res[j-1])

            pl.plot(np.arange(1, max_Cdet+1, 1), res, linewidth=3,
                    color = col[k], linestyle='-',
                    label='HEMT C = {:.2} F'.format(v.Chemt))
            i += 1
            k += 1

        '''
        pl.plot(np.arange(1, max_Cdet+1, 1), res, label='Res eV with C= {:2}
                + 'pF\n'
                .format(v.Chemt*1e12)+'$e_0$={0:.2e} $e_a$={1:.2e}'
                .format(v.e0, v.ea)+' $i_0$={0:.2e} $i_a$={1:.2e}'
                .format(v.i0, v.ia), linewidth=2, color=col[i])
        '''
        if v == Fet_dimitri:
            # modele de dimtri du coup aucun calcul utilisable
            df = 1
            freq = np.arange(1, 50000, 1)
            max_Cdet = 150
            Cfil = 70e-12
            Cload = 10e-12
            res = np.zeros(max_Cdet)

            for Cd in np.arange(1, max_Cdet+1, 1):
                # on calcul tout ici
                j = Cd
                edac = 2.32e-8
                Z = (mbh.Z_c(freq, Cd * 1e-12) ** (-1) +
                     mbh.Z_c(freq, Cfil) ** (-1) +
                     mbh.Z_c(freq, Cload) ** (-1)) ** (-1)
                noise_f = np.sqrt(v.en_(freq) ** 2 +
                                  abs((v.in_(freq) +
                                      ((edac / mbh.Z_c(freq, Cload)))
                                       ) * Z) ** 2)

                res[j-1] = mbh.resolution_f(noise_f**2, abs(Z))

            pl.plot(np.arange(1, max_Cdet+1, 1), res, linewidth=3,
                    linestyle='-', color='gold',
                    label='FET'.format(v.Chemt*1e12))
            i += 1

    pl.xlabel('Cdetector[pF]', fontsize=24)
    pl.ylabel('Resolution [eV]', fontsize=24)
    pl.xticks(fontsize=22)
    pl.yticks(fontsize=22)
    # pl.axis([0, 100, 0, 50])
    NiceGrid()
    pl.legend(loc='upper right', fontsize=18)
    figManager = pl.get_current_fig_manager()
    figManager.window.showMaximized()
    pl.show()


def figure_compa(hemt):
    """
    Figure montrant toutes les contributions pour le bruit du Hemt

    """
    pl.figure('compa')

    pl.title('Impedance', fontweight='bold')
    pl.suptitle('Tb={0:}K Rb={1:.2e}$\\Omega$  Cd={2:}F'
                .format(*compo)
                + ' Cp={3:}F Cc={4:}F Cfb={5:}F'
                .format(*compo),
                y=0.92, fontsize=10)

    plot = 1

    for v in hemt:

        # Calcul du bruit total
        ax = pl.subplot(2, 2, plot)
        print((v.Chemt))
        noise = mbh.total_noise(freq, v, Tb, Rb, Cd, Cp, Cc, Cfb, ifb)
        i = 0

        # Plot du bruit total
        pl.loglog(freq, noise, label='bruit total', linewidth=2, color='black')
        i += 1
        # Calcul pour les dif contri
        Zb = (mbh.Z_b(freq, Rb, Cd, Cp, Cc, Cfb, v.Chemt))

        contri_in = v.in_(freq)*np.abs(mbh.Z_n(freq, Rb,
                                               Cd, Cp, Cc, Cfb, v.Chemt))

        contri_en = v.en_(freq)

        res_f = mbh.resolution_f(noise**2, np.abs(Zb))

        # Plot du signal en frequentiel

        contri_ib = np.abs(Zb) * mbh.ejohnson(Tb, Rb) / Rb

        pl.loglog(freq, contri_in, color='black', linestyle=':',
                  label='contribution in')

        pl.loglog(freq, contri_en, color='black', linestyle='--',
                  label='contribution en \n res={:.3}eV '.format(res_f))
        pl.loglog(freq, contri_ib, color=col[i], label='contribution de Rb')

        pl.xlabel('Frequency [Hz]')

        pl.ylabel('Noise LPSD [V/$\\sqrt{Hz}$]')

        NiceGrid()
        text = ('{:2}pF Hemt'.format(v.Chemt*1e12))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)

        # place a text box in upper left in axes coords
        ax.text(0.1e4, 5e-8, text, fontsize=10,
                verticalalignment='top', bbox=props)
        pl.axis([0.9, 2e4, 1e-10, 1e-7])
        pl.legend(loc='lower left')
        plot += 1

    figManager = pl.get_current_fig_manager()
    figManager.window.showMaximized()
    pl.show()


# liste des hemts avec leurs cara propres
hemt200 = mbh.HEMT(0.18, 5.2, 0, 4.5e-5, 21, 0, 236)
hemt100 = mbh.HEMT(0.23, 6.7, 0, 1.7e-4, 15.5, 0, 103)
hemt40 = mbh.HEMT(0.12, 14.9, 0, 5.0, 7.59, 0, 33)
hemt4 = mbh.HEMT(0.21, 37, 0, 4.3e-6, 2.21, 0, 4.6)
hemt2 = mbh.HEMT(0.4, 91.4, 0, 3.1, 1.8, 0, 1.8)
Fet_dimitri = mbh.HEMT(1.61, 9.61, 30.8, 78.7, 48.6, 0.311, 15)
hemt_alex100 = mbh.HEMT(0.22, 7.3, 0, 0, 16, 0, 100)
hemt_alex36 = mbh.HEMT(0.12, 16.6, 0, 0, 9, 0, 36)
hemt_alex4 = mbh.HEMT(0.21, 44, 0, 0, 2.2, 0, 4.6)

hemt_list = [hemt200, hemt100, hemt40, hemt4]
# hemt_list = [hemt100]


[Tb, Rb, Cd, Cp, Cc, Cfb] = [20e-3, 1e10,
                                         20e-12, 1e-12, 2e-9, 1e-12]
compo = mbh.composant(20e-3, 1e10, 20e-12, 1e-12, 2e-9, 1e-12)
ifb = 0

fmax = 50000

freq = np.arange(1, fmax+1, 1)

# des courbes
# figure([hemt4])

# figure_impedance(hemt_list)
# figure_compa(hemt_list)

hemt_list = [hemt200, hemt100, hemt40, hemt_alex4, hemt2, Fet_dimitri]
fig_scan_Cdet_compa(hemt_list)
