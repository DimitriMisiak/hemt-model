#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:11:36 2019

@author: filippini

TEST UPROOT
"""

import numpy as np 
import uproot
import matplotlib.pyplot as pl
import scipy.signal as sgl



def PSD_unroot():

    fp = ('/home/filippini/Documents/DATA/tc04l004/RED50/' +
          'ProcessedData_tc04l004_S00_RED50_ChanTrig0.root')

    A = uproot.open(fp)

    B = A['RunTree_Normal']

    PSD_filt = B['PSD_Filt'].array()

    filter_order = B['filter_order'].array()

    cutoff_freq = B['cutoff_freq'].array()

    n = filter_order[0]

    fc = cutoff_freq[0]

    f = np.arange(0, 201)

    Gain = (1 + (fc / f) ** (2 * n)) ** (-0.5)

    PSD_filt = PSD_filt.T[0][0:201].T[0] / Gain

    Gain_adu = B['Chan_Gain'].array()[0][0]
    PSD_filt = PSD_filt * Gain_adu * 1e-9
    #
    #
    # test calcul resolution
    pulse = B['Template_Heat_Raw'].array()
    pulse = pulse[0].T[0][0:201]
    pulse = pulse 
    nep2 = pulse[1:201] ** 2 /PSD_filt[1:] ** 2
    
    of = B['OF_Filt'].array()
    of = of[0].T[0][0:201] * Gain
    # pl.loglog(of)
    res = np.sum(nep2[:])
    res = (res ** (-0.5)) / np.sqrt(2)
    print(B['OF_Filt_resolutions'].array()[0][0], 'verif calcul', res)
    #
    #
    #
    return PSD_filt  , pulse  ** 2
