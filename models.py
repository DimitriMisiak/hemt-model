#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Different models used to fit the thermal pulses.

Author:
    Dimitri Misiak (misiak@ipnl.in2p3.fr)
"""


import numpy as np


def physicond(amp, t_phi, t_th):
    """ Function enforcing the physical conditions of the system. Namely,
    the perturbation in temperature cannot be a negative one. Thus function
    assures a positive or null amplitude,
    positive or null characteristic times,
    and quicker rising exponential than the decaying exponential.

    Parameters
    ----------
    amp : float
        Amplitude of the exponential pulse.
    t_phi : float
        Characteristic time of the decaying exponential.
    t_th : float
        Characteristic time of the rising (thermalisation) exponential.

    Returns
    -------
    Same as parameters with physical constraints applied.
    """
    # no negative amplitude
    if amp <= 0:
        amp = 0

    # no negative char. time
    if t_th <= 0:
        t_th = 1e-20

    # rising time shorter than decay time
    if t_phi <= t_th:
        t_phi = t_th

    return amp, t_phi, t_th


def model_1exp(a, t, s, t0=0):
    """Returns the decaying exponential function with heaviside
    and thermalisation exponential which characteristics
    are given by the arguments.

    Parameters
    ----------
    a : float
        Amplitude of the exponential.
    t : float
        Characteristic time of the exponential.
    s : float
        Thermalisation time of the phonons.
    t0 : float, optional
        Time offset of the function.

    Returns
    -------
    aux : function
        Takes as arguments an array_like for the time in s.
    """
    a, t, s = physicond(a, t, s)

    # auxilary function taking as argument a time array.
    def aux(array_t, details=False):

        # applying time offset
        t_var = array_t - t0

        # model expression coming from symbolic calculation
        pulse = np.heaviside(t_var, 1.)
        ind = np.where(pulse != 0.)
        pulse[ind] = a * (np.exp(-t_var[ind]/t)-np.exp(-t_var[ind]/s))

        # same behavior as model_2exp and model_3exp
        if details is True:
            return pulse, pulse
        elif details is False:
            return pulse

    return aux


def model_1exp_fft(ns, a, t, s, t0=0):
    """ FREQUENCY SPACE !
    Returns the decaying exponential function with heaviside
    and thermalisation exponential which characteristics
    are given by the arguments.

    Parameters
    ----------
    ns : float
        FFT Normalization factor. Should be the number of point
        in time space.
    a : float
        Amplitude of the exponential.
    t : float
        Characteristic time of the exponential.
    t0 : float, optional
        Time offset of the function.
    thermal : boolean, optional
        Add the thermalisation exponential term to the function.
    s : float
        Thermalisation time of the phonons.

    Returns
    -------
    aux : function
        Takes as arguments an array_like for the time in s.
    """
    a, t, s = physicond(a, t, s)

    # auxilary function taking as argument a time array.
    def aux(nu_array):
        rising = a*s / (1j * 2*np.pi * s*nu_array + 1)
        decaying = a*t / (1j * 2*np.pi * t*nu_array + 1)
        offset_phase = np.exp(-1j * 2*np.pi * t0*nu_array)
        return (decaying - rising) * offset_phase * ns

    return aux


def model_2exp(a1, a2, t1, t2, s, **kwargs):
    """ Returns a function which is the linear combination of 3 decaying
    exponentials which characteristics are given by the arguments.

    Parameters
    ----------
    (a1, a2, a3) : array_like of floats
        Amplitudes of the exponentials.
    (t1, t2, t3) : array_like of floats
        Characteristic time of the exponentials.
    Additional kwargs will be passed to the model_1exp function. Refer to it
    to check for the keywords.

    Returns
    -------
    aux : function
        Takes as arguments an array_like for the time in s. The argument
        details modifies the return of aux so that it gives each exponential
        separated.

    See also
    --------
    model_1exp
    """
    # auxilary function taking a time array as argument.
    def aux(array_t, details=False):
        exp1 = model_1exp(a1, t1, s, **kwargs)(array_t)
        exp2 = model_1exp(a2, t2, s, **kwargs)(array_t)
        exp_tot = exp1+exp2

        # can return exponential components
        if details is True:
            return exp_tot, exp1, exp2
        elif details is False:
            return exp_tot

    return aux


def model_2exp_fft(ns, a1, a2, t1, t2, s, **kwargs):
    """ FREQUENCY SPACE
    Returns a function which is the linear combination of 3 decaying
    exponentials which characteristics are given by the arguments.

    Parameters
    ----------
    ns : float
        FFT Normalization factor. Should be the number of point
        in time space.
    (a1, a2, a3) : array_like of floats
        Amplitudes of the exponentials.
    (t1, t2, t3) : array_like of floats
        Characteristic time of the exponentials.
    Additional kwargs will be passed to the model_1exp function. Refer to it
    to check for the keywords.

    Returns
    -------
    aux : function
        Takes as arguments an array_like for the time in s. The argument
        details modifies the return of aux so that it gives each exponential
        separated.

    See also
    --------
    model_1exp
    """
    # auxilary function taking a time array as argument.
    def aux(nu_array, details=False):
        exp1 = model_1exp_fft(ns, a1, t1, s, **kwargs)(nu_array)
        exp2 = model_1exp_fft(ns, a2, t2, s, **kwargs)(nu_array)
        exp_tot = exp1 + exp2

        # can return exponential components
        if details is True:
            return exp_tot, exp1, exp2
        elif details is False:
            return exp_tot

    return aux


def model_3exp(a1, a2, a3, t1, t2, t3, s, **kwargs):
    """Returns a function which is the linear combination of 3 decaying
    exponentials which characteristics are given by the arguments.

    Parameters
    ----------
    (a1, a2, a3) : array_like of floats
        Amplitudes of the exponentials.
    (t1, t2, t3) : array_like of floats
        Characteristic time of the exponentials.
    Additional kwargs will be passed to the model_1exp function. Refer to it
    to check for the keywords.

    Returns
    -------
    aux : function
        Takes as arguments an array_like for the time in s. The argument
        details modifies the return of aux so that it gives each exponential
        separated.

    See also
    --------
    model_1exp
    """
    # auxilary function taking a time array as argument.
    def aux(array_t, details=False):
        exp1 = model_1exp(a1, t1, s, **kwargs)(array_t)
        exp2 = model_1exp(a2, t2, s, **kwargs)(array_t)
        exp3 = model_1exp(a3, t3, s, **kwargs)(array_t)
        exp_tot = exp1+exp2+exp3

        # can return exponential components
        if details is True:
            return exp_tot, exp1, exp2, exp3
        elif details is False:
            return exp_tot

    return aux


def model_3exp_fft(ns, a1, a2, a3, t1, t2, t3, s, **kwargs):
    """Returns a function which is the linear combination of 3 decaying
    exponentials which characteristics are given by the arguments.

    Parameters
    ----------
    ns : float
        FFT Normalization factor. Should be the number of point
        in time space.
    (a1, a2, a3) : array_like of floats
        Amplitudes of the exponentials.
    (t1, t2, t3) : array_like of floats
        Characteristic time of the exponentials.
    Additional kwargs will be passed to the model_1exp function. Refer to it
    to check for the keywords.

    Returns
    -------
    aux : function
        Takes as arguments an array_like for the time in s. The argument
        details modifies the return of aux so that it gives each exponential
        separated.

    See also
    --------
    model_1exp
    """
    # auxilary function taking a time array as argument.
    def aux(nu_array, details=False):
        exp1 = model_1exp_fft(ns, a1, t1, s, **kwargs)(nu_array)
        exp2 = model_1exp_fft(ns, a2, t2, s, **kwargs)(nu_array)
        exp3 = model_1exp_fft(ns, a3, t3, s, **kwargs)(nu_array)
        exp_tot = exp1 + exp2 + exp3

        # can return exponential components
        if details is True:
            return exp_tot, exp1, exp2, exp3
        elif details is False:
            return exp_tot

    return aux


def dim_ifft(fft_funk, fs, ns):
    """ Numpy friendly function meant to convert the fft models
    into time space (temporal pulse).

    Parameters
    ----------
    fft_funk : function
        Fourier transform model function.
    fs : int
        Sampling frequency.
    ns : int
        Numbers of points in time space (should be 2*frequencies).

    Returns
    -------
    time_pulse : array
        Array containing the pulse in time space.
    """
    f_array = np.fft.fftfreq(ns, fs**-1)
    time_pulse = np.fft.ifft(fft_funk(f_array))

    return time_pulse


def dim_psd(dft, timestep, timelength):
    """ Computes the Power Spectral Density from
    the Discrete Fourier Transform.

    Parameters
    ----------
    dft : array
        DFT array.
    timestep : float
        Time step.
    timelength : float
        Time length.

    Returns
    -------
    psd : array
        PSD array in A^2/Hz, where A is the unit of the signal.
    """
    psd_corr = 2 * timestep**2 / timelength
    psd = np.abs(dft)**2 * psd_corr

    return psd
