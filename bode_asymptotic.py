'''This programme calculates and plots the asymptotic log-magnitude response 
of LTI systems without poles or zeros at the origin

Published under GNU GPL v3
Copyright (c) 2016 Tanmoy Dasgupta (thetdg@live.com)

This is a free software (free, as in freedom).
However, there is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.

You can share/modify this programme 
as long as you keep this copyright notice intact
'''


from __future__ import division

from control.matlab import bode, tf
from numpy import (poly, logspace, log10, amax, amin, roots, sort, floor, 
                    zeros, inf, zeros_like, cumsum)
from matplotlib.pyplot import (close, subplots, semilogx, title, 
                    xlabel, ylabel, axis, grid, show)
from matplotlib.ticker import AutoMinorLocator


def plot_range(num, den):

    # The corner frequencies
    zero = sort(abs(roots(num)))    
    pole = sort(abs(roots(den)))
    
    # Calculate the minimum and maximum corner frequencies needed
    if len(pole) == 0:
        corner_min = zero[0]
        corner_max = zero[-1]
    
    elif len(zero) == 0:
        corner_min = pole[0]
        corner_max = pole[-1]
        
    elif len(zero) > 0 and len(pole) > 0:
        corner_min = min(zero[0], pole[0])
        corner_max = max(zero[-1], pole[-1]) 
    
    else:
        corner_min, corner_max = 0.1, 10
    
    # start from 2 decades lower than the lowest corner 
    # and end at 2 decades above the highest corner
    freq_range = [10 ** (floor(log10(corner_min)) - 1), 
                   10 ** (floor(log10(corner_max)) + 2)]
    
    return freq_range


def asymptote(num, den):

    # create a Python list for the zeros and the poles of the system
    zero = list(sort(abs(roots(num))))    
    pole = list(sort(abs(roots(den))))
    
    #calculate the low frequency gain -- type 0 system
    lf_gain = 20 * log10(abs((num[-1] / den[-1])))
    
    # create an empty matrix to contain the corner frequencies and
    # the corresponding slope indicator (+1 or -1)
    corners = zeros((len(zero) + len(pole) + 2, 2))
    
    starting_freq, end_freq = plot_range(num, den)
    corners[0] = [starting_freq, 0]
    corners[-1] = [end_freq, 0]
    
    # take the first elements from the list of poles and zeros
    # compare them and assign the slope indicator
    # delete the corresponding valuefrom the original list of poles and zeros
    for count in range(len(zero) + len(pole)): 
            
        if len(zero) > 0:
            a = zero[0]
        else:
            a = inf
            
        if len(pole) > 0:
            b = pole[0]
        else:
            b = inf
        
        c = min(a, b)
                
        if c == a:
            corners[count + 1] = [c, 1]
            if len(zero) > 0:
                zero.pop(0)

        if c == b:
            corners[count + 1] = [c, -1]
            if len(pole) > 0:
                pole.pop(0)
         
    # now calculate the gains at the corners using 
    # gain = +/- 20log10(upper_corner / lower_corner) 
    asymptotic_gain = zeros_like(corners)   
    asymptotic_gain[0, 1] = lf_gain
    asymptotic_gain[1, 1] = lf_gain
    
    gain = lf_gain
    multiplier = cumsum(corners[:, 1])
    for k in range(2, len(corners)):    
        gain += multiplier[k-1] * 20 * log10(corners[k, 0] / corners[k-1, 0])
        asymptotic_gain[k, 1] = gain
        
    asymptotic_gain[:, 0] = corners[:, 0]
                   
    return asymptotic_gain


def plot_config(mag, omega, mag_as, omega_as):
    
    # plot configurations
    # change it as you need
    close('all')
    pad = 1
    fig, ax = subplots(figsize=(15, 7))
    
    # These are the plots
    semilogx(omega, mag, lw=2)
    semilogx(omega_as, mag_as, 'g', lw=2)
    
    # Blah blah
    hfont = {'fontname':'cmr10'}
    title(r'Asymptotic log-magnitude Bode plot', fontsize=24, **hfont)
    xlabel('Frequency $\omega$ in rad/s', fontsize=24, **hfont)
    ylabel('Magnitude $|G(j\omega)|$ in dB', fontsize=24, **hfont)
    axis([amin(omega), amax(omega), 
            round(min(amin(mag), amin(mag_as)) - pad), 
            round(max(amax(mag), amax(mag_as)) + pad)])
    minor_locator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minor_locator)
    ax.tick_params(axis='both', labelsize=24)
    grid(b=True, which='major', color='k', linestyle='-')
    grid(b=True, which='minor', color='r', alpha=.5, linestyle='-')
    ax.set_xticklabels(ax.get_xticks(), hfont)
    ax.set_yticklabels(ax.get_yticks(), hfont)
    show()


def bode_as(num, den, freq_range=None):
    '''Creates the asymptotic log-magnitude Bode plot for type zero transfer 
    functions.
    
    Example : for the system G(s) = (4s + 1) / (3s + 2), use
    bode_as([4, 1], [3, 2])
    
    for the system G(s) = [(s + 1)(s + 20)] / [(s + 2)(s + 8)], use
    bode_as(poly((-1, -20)), poly((-2, -8)))
    '''
    
    G = tf(num, den)
    
    # If the user provides the range of frequencies, take it
    # Else calculate the most suitable range
    if freq_range:
        freq_min, freq_max = freq_range
    else:
        freq_min, freq_max = plot_range(num, den)
        
    omega = logspace(log10(freq_min), log10(freq_max), 1000)
    
    # The actual Bode plot calculation
    mag, phase, freq = bode(G, omega, dB=True, Plot=False)
    
    # Asymptotic log-magnitude Bode plot calculation
    asymptotic_plot_data = asymptote(num, den)
    omega_as, mag_as = asymptotic_plot_data[:, 0], asymptotic_plot_data[:, 1]
    
    plot_config(mag, omega, mag_as, omega_as)
