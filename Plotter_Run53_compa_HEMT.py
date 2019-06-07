#!/usr/bin/python

# coding:utf-8


####################################

##

#  Compa hemt modele et run en voie chaleur

##

####################################

import numpy as np

import matplotlib.pyplot as pl

import os

import time

import copy

import logging

import sys

import re

import glob

import Modele_ampli_Hemt_bruit as mbh

import unroot_test as root

# impedance du syteme avant le HEMT

def Zc(freq, C):
    return 1/(2*1j*np.pi*freq*C)


def Z(freq, R, C):
    return abs(1/((1/R)+(1/Zc(freq, C))))


def Zfet(freq, R, Cfil, Cload):
    return (R ** (-1) + Zc(freq, Cfil) ** (-1) + Zc(freq, Cload) ** (-1)) ** (-1)


def ejohnsonfet(T, R):
    return (4*kb*T*R)


def ejohnson(freq, T, R, C):
    return abs((Zc(freq, C)/(Zc(freq, C)+R))*(4*kb*T*R))


composant = [Tb, Rb, Cd, Cp, Cc, Cfb] = [20e-3, 1e10,
                                         2.5e-11, 1e-11, 2e-9, 1e-12]

t200 = mbh.HEMT(0.18, 5.2, 0, 4.5e-5, 21, 0, 236)

t100 = mbh.HEMT(0.23, 6.2, 0, 1.7e-4, 15.5, 0, 103)

t40 = mbh.HEMT(0.12, 14.9, 0, 5.0, 7.59, 0, 33)

t4 = mbh.HEMT(0.21, 37.1, 0, 4.3e-6, 2.21, 0, 4.6)

t2 = mbh.HEMT(0.4, 91.4, 0, 3.1, 1.8, 0, 1.8)

Fet_dimitri = mbh.HEMT(1.61, 9.61, 30.8, 78.7, 48.6, .311, 70)

Fet_alex = mbh.HEMT(0.5, 7.3, 0, 20, 50, 0, 50)

Hemt_list = [t2, t4, t40, t100, t200]

kb = 1.38e-23

save = len(sys.argv)

t0 = time.time()

# ##palette couleur ###

col = ['black', 'orangered', 'royalblue', 'fuchsia', 'orange', 'darkcyan']

col10 = pl.cm.Set1(np.linspace(0, 1, 9))


axiss = [1, 200, 0.5e-9, 100e-9]  # axes pour tout les plots

# ### MATRICE D'INFO DES RUNS, COMMENTER LES RUNS A NE PAS PLOTTER ###


run22 = [['rc27g000', '40 mK', '40 kOhms', '45 kOhms', 'C0', '-'],
         ['rc27g001', '35 mK', '90 kOhms', '110 kOhms', 'C1', '-'],
         ['rc27g002', '30 mK', '250 kOhms', '300 kOhms', 'C2', '-'],
         ['rc27g003', '25 mK', '865 kOhms', '920 kOhms', 'C3', '-'],
         ['rc27g005', '20 mK', '3.7 MOhms', '3.3 MOhms', 'C4', '-'],
         ['rc28g001', '18 mK', '8.3 MOhms', '6.8 MOhms', 'C5', '-'],
         ['rc28g002', '16 mK', '15.6 MOhms', '11 MOhms', 'C6', '-'],
         ['rc28g003', '14 mK', '28.6 MOhms', '15.6 MOhms', 'C7', '-'],
         ['rc28g004', '12 mK', '45.6 MOhms', '21.8 MOhms', 'C8', '-']]


MY_FORMAT = "%(asctime)s.%(msecs)03d %(levelname)-6s %(message)s %(funcName)s"

logging.basicConfig(format=MY_FORMAT, datefmt='%H:%M:%S')

MY_LOGGER = logging.getLogger()

MY_LOGGER.setLevel(logging.INFO)


datapath = '/home/filippini/Documents/Run53/DATA'  # path des datas

datapath2 = '/home/filippini/Documents/Run53'  # path du run


# ## Extraction des noms des fichiers et des variables ####

path_npy = glob.glob(datapath2+'/npy/*.npy')
# cherche tout les fichiers npy pour extraire les infos utile

path_npy = np.sort(path_npy)

test = np.load(path_npy[0])  # charge les donnees du npy dans un tableau

size = len(test)

j = 0

filename = copy.copy(path_npy)

for i in path_npy:

    filename[j] = os.path.basename(i)

    filename[j] = os.path.splitext(filename[j])[0]
# Nom du fichier contenant le run la frequence d echan et la temperature
    j += 1

# initialisation des variables
T = list(range(j))

fs = list(range(j))

run = list(range(j))

sub = list(range(j))

R = np.zeros((size, j))

name = list()


name = re.findall("Hz_([0-9a-zA-Z]+)", filename[0])

# recupere le nom des detecteurs
for i in range(size-1):

    name.append(re.findall("MO_([0-9a-zA-Z]+)", filename[i]))
    name[i+1] = str(name[i+1]).strip("[]'")

# recupere toutes les donnees utile
for i in range(j):

    T[i] = re.findall("([.0-9]+)mK", filename[i])  # temperature

    fs[i] = re.findall("fs=([0-9]+)", filename[i])  # frequence d echanti

    run[i] = re.findall("(^[a-z0-9]+)", filename[i])  # nom du run

    sub = re.findall("([0-9.]+)MO", filename[i])  # resistance
    # Pour chaque detecteur recupere la resistance du ntd
    for a in range(size):

        sub[a] = str(sub[a]).strip("[]")

        R[a][i] = sub[a]

    T[i] = str(T[i]).strip("[]'")

    fs[i] = str(fs[i]).strip("[]'")

    run[i] = str(run[i]).strip("[]'")


T = [float(i) for i in T]

fs = [int(i) for i in fs]

# R = [float(i) for i in R]

f1 = range(251)

f2 = np.arange(1, 50000, 1)

pi = 3.14

save = 100

MY_LOGGER.info('init script fait, figs ongoing')


def NicePlot():

    pl.xlabel('Freq [Hz]', fontsize=24)
    pl.ylabel('Noise [V/$\\sqrt{Hz}$]', fontsize=24)
    pl.grid(b=True, which='major', color='black', linestyle='-')
    pl.grid(b=True, which='minor', color='silver', linestyle=':')


def NiceGrid():

    pl.grid(b=True, which='major', color='black', linestyle='-')

    pl.grid(b=True, which='minor', color='silver', linestyle=':')


nom = 0


def Fig32():

    path_txt = glob.glob(datapath2+'/DATA/Run32/*.txt')

    path_txt = np.sort(path_txt)

    freq, PSDmatrix = np.loadtxt(path_txt[0], unpack=True)
    pl.loglog(freq, PSDmatrix, color=col10[0], linestyle='-')


def Fig2(nom):  #

    # boucle sur le nom des detecteurs
    for i1 in [0]:

        Fig = pl.figure('Fig'+str(nom), figsize=(12.8 / 1.2, 8 / 1.2))

        nom += 1

        Fig.suptitle('FET experimental '+str(name[i1])+'vs HEMT model Noise',
                     fontsize=12, fontweight='bold')

        j = 0
        j1 = 0

        path_txt = glob.glob('/home/filippini/Documents/Run53/DATA/Run32/*.txt')

        path_txt = np.sort(path_txt)

        freq, PSDmatrix32 = np.loadtxt(path_txt[4], unpack=True)

        # boucle sur  la frequence d echantillonnage
        for k in [1]:

            j = k

            ax = pl.subplot(1, 1, 1)

            NiceGrid()

            pl.ylabel('Noise [V/$\\sqrt{Hz}$]', fontsize=24)

            pl.xlabel('Frequency [Hz]', fontsize=22)

            a = 4
            a += k

            # boucle sur les temperatures
            for i in [a]:
                labela = str(run[i])+'_fs='+str(fs[i])+'Hz_T='+str(T[i])+'mk'

                PSDmatrix = np.load(path_npy[i])
                '''
                pl.loglog(range(int(np.array(fs[i]/2+1))),
                          np.sqrt(PSDmatrix[k, :]), color=col10[j1],
                          linestyle='-', label=labela, linewidth=2)
                                
                pl.text(0.15+0.24*j, 0.05, str(T[i]) + 'mK ' + str(R[0][i])
                        + 'MO', horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes, color=col10[j1],
                        fontsize=8, fontweight='bold')
                '''
                j1 += 1

            j += 2

        j -= 2

        en = np.zeros(int(np.array(fs[a])))
        entot = np.zeros(int(np.array(fs[a])))
        tab_in = np.zeros(int(np.array(fs[a])))
        impedance = np.zeros(int(np.array(fs[a])))
        bjohnson = np.zeros(int(np.array(fs[a])))
        idac = np.zeros(int(np.array(fs[a])))
        tes = 0

        i2 = 0
    # calcul hemt
    for v in [t100]:
        for l in range(1, int(np.array(fs[a]))):

            impedance[l-1] = Z(l, R[0][i2]*1e6, v.Chemt)
            en[l-1] = v.en_(l)

            tab_in[l-1] = v.in_(l)

            bjohnson[l-1] = ejohnson(l, T[i2]*1e-3, R[0][i2]*1e6, v.Chemt)

            entot[l-1] = np.sqrt(en[l-1]**2 +
                                 (tab_in[l-1]*impedance[l-1])**2+bjohnson[l-1])

            tab_in[l-1] *= impedance[l-1]

        pl.loglog(range(int(np.array(fs[a]))), entot, linestyle='-',
                  color=col[tes], linewidth=4, label='Bruit total HEMT')
#        pl.loglog(freq, PSDmatrix32, color=col10[2], linestyle='-',
#                  label='RUN 22')

        pl.loglog(range(int(np.array(fs[a]))), tab_in, linestyle=':',
                  color=col[tes], label='Bruit en courant')

        pl.loglog(range(int(np.array(fs[a]))), en, linestyle='--',
                  color=col[tes], label='Bruit en tension')

        pl.loglog(range(int(np.array(fs[a]))), np.sqrt(bjohnson),
                  linestyle='-.', color=col[tes],
                  label='Bruit détecteur')
        pl.axis(axiss)
        pl.xticks(fontsize=22)
        pl.yticks(fontsize=22)
        pl.legend()
        pl.legend(loc='upper right', fontsize=18)
       
        
        
        
        # calcul resolution
        noise = mbh.total_noise(freq[:201], v, T[i2]*1e-3, R[0][i2]*1e6,
                                Cd, Cp, Cc, Cfb)
        pulse = mbh.pulse_t(fs[a])
        PSD04, PSD_signal = root.PSD_unroot()
        res = (mbh.resolution_t(PSDmatrix[k,:], PSD_signal, fs=fs[a])
              * 5890 / 644)
        print("resolution PSD {:2}".format(res))
#        res = mbh.resolution_t(entot**2, pulse, fs=fs[a], i_range=200)
#        print("resolution hemt {:2}".format(res))
        res = mbh.resolution_t(noise**2, PSD_signal, fs=fs[a])/np.sqrt(2) * 5890 / 644
        print("resolution hemt {:2} new formule".format(res))
        res04 = (mbh.resolution_t(PSD04, PSD_signal, fs=400)/np.sqrt(2))* 5890 / 644
        print("resolution PSD04 {:2}".format(res04))
        tes += 1





    Rfet = R[0][a]
    Ta = T[i2] * 1e-3
    v = Fet_dimitri
    edac = 2.32e-8
    Cfil = 1e-11
    # Calcul fet
    for l in range(1, int(np.array(fs[a]))):

        impedance[l-1] = np.abs(Zfet(l, Rfet*1e6, Cfil, v.Chemt))

        en[l-1] = v.en_(l)

        tab_in[l-1] = v.in_(l)

        ijohnson = np.sqrt(ejohnsonfet(Ta*1e-3, Rfet*1e6)) / (Rfet * 1e6)

        idac[l-1] = edac / np.abs(Zc(l, 10e-12))

        entot[l-1] = np.sqrt(en[l-1] ** 2 +
                             (tab_in[l-1] ** 2 + idac[l-1] ** 2 + ijohnson ** 2) *
                             (impedance[l-1]) ** 2)

        tab_in[l-1] *= impedance[l-1]
        # i2+=2
#    res = mbh.resolution_t(entot**2, pulse, i_range=400) *  5890 /644 
#    print("resolution fet {:2}".format(res))
    '''
    pl.loglog(range(int(np.array(fs[a]))), entot, linewidth=5, linestyle=':',
              color=col[tes+3],
              label='Fet pour R='+str(Rfet)+'MO C='+str(v.Chemt*1e12)+'pF')
    pl.loglog(range(int(np.array(fs[a]))), tab_in, linestyle=':',
              color='red',
              label='in pour R='+str(Rfet)+'MO C='+str(v.Chemt*1e12)+'pF')

    pl.loglog(range(int(np.array(fs[a]))), ijohnson * impedance, linestyle=':',
              color=col[tes],
              label='bruit johnson')

    pl.loglog(range(int(np.array(fs[a]))), en, linestyle='--',
              color=col[tes], label='en pour R=' +
                                    str(R[0][i2]) + 'MO C=' +
                                    str(v.Chemt*1e12) + 'pF')
    '''
#    pl.axis(axiss)
#    pl.legend()
#    pl.legend(loc='upper right', fontsize='x-small')

    j = 0

    if save == 2:
        pl.savefig(datapath2 + '/' +
                   'RUN53_Bilan_evolution_temp'+str(i1+1)+'.pdf')

        pl.savefig(datapath2 + '/' +
                   'RUN53_Bilan_evolution_temp'+str(i1+1)+'.png')

    MY_LOGGER.info('DONE')

    return nom


save = 0

pl.close("all")

# Fig32()
datapath2 += '/plot'
nom = Fig2(nom)

print('temps script', time.time() - t0)

pl.show()
