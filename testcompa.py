#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 16:51:04 2019

@author: filippini

Comparaison des spectres de bruits hemts
"""


import numpy as np

import matplotlib.pyplot as pl

import Modele_ampli_Hemt_bruit as mbh


def nicegrid():
    """
    Grille qui est en pointille en loglog
    """
    pl.grid(b=True, which='major', color='black', linestyle='-')

    pl.grid(b=True, which='minor', color='silver', linestyle=':')


t100 = mbh.HEMT(0.23, 6.2, 0, 1.7e-4, 15.5, 0, 103)
DATA_PATH = '/home/filippini/Documents/DATA/RUN55/txt/'
INFO_RUN = 'Rien'

# point 100mV 1m1 voie B
#PATH_TXT = [['PSD_matrix_20190509_17h56.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 150mV voie A'],
#            ['PSD_matrix_20190509_17h53.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 150mV voie B']
#            ]

#PATH_TXT = [['PSD_matrix_20190510_12h21.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 100mV voie A'],
#            ['PSD_ma'trix_20190510_12h26.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 100mV voie B']
#            ]

#PATH_TXT = [['PSD_matrix_20190510_13h02.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 200mV voie A'],
#            ['PSD_matrix_20190510_13h05.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 200mV voie B']
#            ]

#PATH_TXT = [['PSD_matrix_20190510_15h56.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 100mV voie A'],
#            ['PSD_matrix_20190510_16h00.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 100mV voie B']
#            ]

#PATH_TXT = [['PSD_matrix_20190510_13h02.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 200mV voie A gain 7.87'],
#            ['PSD_matrix_20190509_17h56.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 150mV voie A gain 6.81'],
#            ['PSD_matrix_20190510_12h21.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 100mV voie A gain 6.37'],
#            ['PSD_matrix_20190510_15h56.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 100mV voie A gain 4.33']
#            ]

#PATH_TXT = [['PSD_matrix_20190510_13h05.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 200mV voie B gain 7.48'],
#            ['PSD_matrix_20190509_17h53.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 150mV voie B gain 6.43'],
#            ['PSD_matrix_20190510_12h26.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '1mA 100mV voie B gain 4.68'],
#            ['PSD_matrix_20190510_16h00.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 100mV voie B gain 4.28']
#            ]

#PATH_TXT = [['PSD_matrix_20190510_16h00.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 100mV voie B gain 4.28'],
#            ['PSD_matrix_20190515_16h54.BIN0_T=4K_fs=100000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 150mV voie B gain 4.44'],
#            ['PSD_matrix_20190515_18h21.BIN0_T=4K_fs=100000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 200mV voie B gain 4.64'],
#            ['PSD_matrix_20190516_10h44.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 75mV voie B gain 3.59']
#            ]    


#PATH_TXT = [['PSD_matrix_20190510_15h56.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 100mV voie A gain 4.33'],
#            ['PSD_matrix_20190515_16h43.BIN0_T=4K_fs=100000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 150mV voie A'],
#            ['PSD_matrix_20190515_18h23.BIN0_T=4K_fs=100000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 200mV voie A gain 4.97'],
#             ['PSD_matrix_20190516_10h41.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 75mV voie A gain 3.85']
#            ]
 
#PATH_TXT = [['PSD_matrix_20190516_15h45.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 150mV 16 mai'],
#            ['PSD_matrix_20190515_16h43.BIN0_T=4K_fs=100000.0Hz_ts=1s_HEMT_Evt',
#            '0.4mA 150mV voie A'],
#            ]

#PATH_TXT = [['PSD_matrix_20190516_16h13.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.7mA 100mV voie B gain 4.87'],
#            ['PSD_matrix_20190516_16h43.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.7mA 150mV voie B gain 5.37'],
#            ['PSD_matrix_20190516_17h04.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.7mA 200mV voie B gain 6.01'],
#             ['PSD_matrix_20190516_17h44.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            '0.7mA 75mV voie B gain 4.36'],
#            ]

PATH_TXT = [
            ['PSD_matrix_20190516_17h41.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
            'Vds=75mV voie A'],
            ['PSD_matrix_20190516_16h10.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
            'Vds=100mV voie A'],
            ['PSD_matrix_20190516_16h40.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
            'Vds=150mV voie A'],
            ['PSD_matrix_20190516_17h01.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
            'Vds=200mV voie A'],
           ]

#PATH_TXT = [['PSD_matrix_20190523_15h45.BIN0_T=4K_fs=200000.0Hz_ts=1s_HEMT_Evt',
#            'relier a rien']]

pl.close('all')
pl.figure("compa", figsize=(12.8 / 1.2, 8 / 1.2))

COL = ['green', 'black', 'orangered', 'purple', 'royalblue', 'fuchsia',
       'orange', 'darkcyan']

pl.title("", fontsize=12)
nicegrid()
for i in range(int(np.size(PATH_TXT)/2)):

    test = np.loadtxt(DATA_PATH+PATH_TXT[i][0])

    pl.loglog(np.sqrt(test), color=COL[i], label=PATH_TXT[i][1],
              linewidth=1)
    
    
freq = np.arange(1,500000,1)
a = t100.en_(freq)
pl.loglog(t100.en_(freq), linewidth=5, linestyle=':', label='Objectif')
pl.ylabel('Noise LPSD[V/$\\sqrt{Hz}$]', fontsize=20)

pl.xlabel('Freq[Hz]', fontsize=20)

pl.legend()

pl.legend(loc='upper right', fontsize=16)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.axis([0.9, 1e4, 1e-10, 1e-7])


def save():
    """
    Super save de la figure
    """
    pl.savefig('/home/filippini/Documents/plot/RUN55/compa'+INFO_RUN+'.png')


# save()
