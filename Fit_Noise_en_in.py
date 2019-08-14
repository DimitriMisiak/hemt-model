#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:24:08 2019

@author: filippini
"""

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

# Read data end extract components with a fit
col = ['black', 'orangered', 'green', 'royalblue', 'fuchsia',
       'orange', 'darkcyan']

def read_data(data_file):
    """
    Extract data from a file txt to an array with the structure
    capacity  
    frequency
    data (frequency)
    """
    
    all_data = np.loadtxt(data_file, delimiter=',')
    
    i = 0
   
    capacity_hemt = []
    
    frequency = []
    
    data = []
    
    # just a little trick to now number of lines
    while i < np.size(all_data.T[:][0]) :
        
        # stock all capacity in a table
        capacity_hemt.append(all_data[i][all_data[i] != 0] * 1e-12)
        
        # stock all frequency in a table
        frequency.append(all_data[i+1][all_data[i+1] != 0])
        
        # stock all data in a table
        data.append(all_data[i+2][all_data[i+2] != 0])

        i += 3
    
    return capacity_hemt, frequency, data

def en_fit(freq, e0, ea):
    """
    Return Noise in tension "en" function of e0 ea and frequency
    """
    return np.sqrt(e0 ** 2 + ea ** 2 / freq + (0 / freq) ** 2)


def in_fit(freq, i0, ia):
    """
    Return Noise in current "in" function of i0 ia and frequency
    """
    return np.sqrt(i0 ** 2 + ia ** 2 * freq + (0 * freq) ** 2)

def fit(frequency, data, component):
    """
    Return parameter of the fit with no errors bars 
    """
    if component == 'en':
        parameter_fit = [curve_fit(en_fit, freq, datas, sigma=0.1 * datas)[0]
                         for datas, freq in zip(data, frequency)]
        
    if component == 'in':
        parameter_fit = [curve_fit(in_fit, freq, datas, sigma=0.1 * datas)[0]
                         for datas, freq in zip(data, frequency)]

    return parameter_fit


def figure (title, capacity, parameters_fit, frequency, data):
    """
    Plot of data with the fit of en and in
    """
    plt.figure(title)
    
    #plot points datas
    for i in np.arange(0, np.size(frequency)):
       
        plt.loglog(frequency[i], data[i], marker = 'o', ls='', 
                   color = col[i])
    
    plt.xlabel('Frequency [Hz]', fontsize=12)

    plt.ylabel('Noise [V/$\\sqrt{Hz}$]', fontsize=12)
    
    # Plot fit en
    if title == 'en':
        
        f_range_full = np.linspace(1, 1000000, 1000000)
        
        for i in np.arange(0, np.size(frequency)):
            plt.loglog(f_range_full, en_fit(f_range_full, *parameters_fit[i]),
                       linestyle='-',
                       label = 'Fit for hemt C = {0:.3e}'.format(capacity[i]) +
                       ' e0 = {0:.3e} ea ='.format(*parameters_fit[i]) +
                       '{1:.3e}'.format(*parameters_fit[i]),
                       color = col[i])
   
    # Plot fit in
    if title == 'in':
        
        f_range_full = np.linspace(1, 1000000, 1000000)
        
        for i in np.arange(0, np.size(frequency)):
            plt.loglog(f_range_full, in_fit(f_range_full, *parameters_fit[i]),
                       linestyle='-',
                       label = 'Fit for hemt C = {0:.3e}'.format(capacity[i]) +
                       ' i0 = {0:.3e} ia ='.format(*parameters_fit[i]) +
                       '{1:.3e}'.format(*parameters_fit[i]),
                       color = col[i])
                       
    plt.grid(b=True, which='major', color='black', linestyle='-')

    plt.grid(b=True, which='minor', color='silver', linestyle=':')

    plt.legend()
# main 

en_file = 'en_fit.txt'
in_file = 'in_fit.txt'


# en

capacity_hemt, en_frequency, en_data = read_data(en_file)

# capacity to plot and have components of fit
capacity_plot = np.array([92, 236, 5]) * 1e-12

for i in np.arange(np.size(en_data)):
    en_data[i] *= 1e-9


# select data

# init
number_table = []

capacity = []

for j in np.arange(0, np.size(capacity_plot)):
    
    for i in np.arange(0, np.size(capacity_hemt)):
        
        if capacity_plot[j] == capacity_hemt[i][0]:
        
            number_table.append(i)
            
            capacity.append(capacity_hemt[i][0])

frequency = []
    
data = []


for i in number_table : 
    
    frequency.append(en_frequency[i])
        
    data.append(en_data[i])

 
parameters_fit = fit(frequency, data, 'en')

figure('en', capacity, parameters_fit, frequency, data)


# in

capacity_hemt, in_frequency, in_data = read_data(in_file)

for i in np.arange(np.size(in_data)):
    in_data[i] *= 1e-18

capacity_plot = np.array([92, 236, 5]) * 1e-12

#select data

number_table = []

capacity = []

for j in np.arange(0, np.size(capacity_plot)):
    
    for i in np.arange(0, np.size(capacity_hemt)):
        
        if capacity_plot[j] == capacity_hemt[i][0]:
        
            number_table.append(i)
            
            capacity.append(capacity_hemt[i][0])

frequency = []
    
data = []


for i in number_table : 
    
    frequency.append(in_frequency[i])
        
    data.append(in_data[i])

 
parameters_fit = fit(frequency, data, 'in')

figure('in', capacity, parameters_fit, frequency, data)
