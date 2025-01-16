#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Erick Jorge Canales-Rodriguez, 2024

from __future__ import division

# Myelin water diffusion models
import numpy as np
from dipy.data import get_fnames
from scipy.special import jv, comb
from dipy.io.gradients import read_bvals_bvecs

#from numpy.math import factorial
from math import factorial
from scipy.special import gamma, gammainc, gammaincc, erf

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
matplotlib.rc('text', usetex = True)

# We assume a trapezoid-shaped diffusion profile
def Radial_signal(Nterms, D, b, radius, BigDelta, smalldelta, E, sin_alpha, Lori, pulse):
    # Create radial signal: assuming the gradient is applied perpendicular
    if pulse == 'rectangle':
        t =  BigDelta - smalldelta/3
    elif pulse == 'trapezoid':
        t = (BigDelta - smalldelta/3) + (E**3)/(30*smalldelta**2) - (E**2)/(6*smalldelta)
    #end
    q      =  np.sqrt(b/t)

    #---------- Lori appproximation
    if Lori==True:
        if pulse == 'rectangle':
            q = q * np.sqrt(t/(BigDelta + smalldelta))
            t = BigDelta + smalldelta
        elif pulse == 'trapezoid':
            q = q * np.sqrt(t/(BigDelta + smalldelta + E))
            t = BigDelta + smalldelta + E
        #end
    #end
    # -----------------------------
    Signal_radial = 0
    for n in range(1, Nterms):
        Signal_radial = Signal_radial + np.exp(-(n**2) * ((D*t)/radius**2) ) *jv(n, radius*q*sin_alpha)**2
    #end
    Signal_radial     = jv(0, radius*q*sin_alpha)**2 + 2*Signal_radial
    return Signal_radial
#end

def SMT_signal(Nterms, D, b, radius, BigDelta, smalldelta, E, Lori, pulse):
    if pulse == 'rectangle':
        t =  BigDelta - smalldelta/3
    elif pulse == 'trapezoid':
        t = (BigDelta - smalldelta/3) + (E**3)/(30*smalldelta**2) - (E**2)/(6*smalldelta)
    #end
    q      =  np.sqrt(b/t)

    #---------- Lori approximation
    if Lori==True:
        if pulse == 'rectangle':
            q = q * np.sqrt(t/(BigDelta + smalldelta))
            t = BigDelta + smalldelta
        elif pulse == 'trapezoid':
            q = q * np.sqrt(t/(BigDelta + smalldelta + E))
            t = BigDelta + smalldelta + E
        #end
    #end
    # -----------------------------

    Signal1 = 0
    Nterms2 = 80
    for k in range(0, Nterms2):
        Ck0 = Coeff(k,0) * (radius*q)**(2*k)
        for j in range(0, k+1):
            Signal1 = Signal1 +  Ck0 * comb(k, j) * Res_Gamma(b, D, j)
    #end

    Signal2 = 0
    for n in range(1, Nterms):
        Exp = np.exp(-(n**2) * (D*t)/(radius**2))
        for k in range(0, Nterms2):
            Ckn = Coeff(k,n) * (radius*q)**(2*(n+k))
            for j in range(0, n+k+1):
                Signal2 = Signal2 + Exp * Ckn * comb(n+k,j) * Res_Gamma(b, D, j)
    #end
    Signal = 0.5*Signal1 + Signal2
    return Signal
#end

def Res_Gamma(b, D, p):
    if b*D == 0: # asymptotic expression for b*D -> 0
        value = ((-1)**p ) * gamma(p + 0.5) / gamma(p + 1.5)
    else:
        value = ((-1)**p ) * ( gamma(p + 0.5) - gamma(p + 0.5)*gammaincc(p + 0.5, b*D) )/(b*D)**(p + 0.5)
    #end
    return value
#end

def Coeff(k,n):
    #value = (-1)**k * (1/(factorial(k)*factorial(2*n+k)))* (1/2)**(2*(n+k)) * factorial(2*(n+k))/(factorial(n+k)*factorial(2*(n+k) - (n+k)))
    value = (-1)**(k) * ( 1/(factorial(k)*factorial(2*n+k)) ) * 0.5**(2*(n+k)) * comb(2*(n+k), n+k)
    return value
#end

def SMT_signal_Gaussian(D, b, radius, BigDelta, smalldelta, E, Lori, pulse):
    t =  BigDelta - smalldelta/3

    #---------- Lori appproximation
    if Lori==True:
        if pulse == 'rectangle':
            t = BigDelta + smalldelta
        elif pulse == 'trapezoid':
            t = BigDelta + smalldelta + E
        #end
    #end
    # ------------------------------------------------------------

    Dr    = (radius**2/(2*t)) * (1 - np.exp(-D*t/(radius**2)))
    value = np.sqrt(np.pi/4) * np.exp(-b*Dr) * erf(np.sqrt(b*(D-Dr))) / np.sqrt(b*(D-Dr))
    return value
#end
