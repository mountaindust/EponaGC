# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 11:56:57 2012

@author: Christopher Strickland
"""

from __future__ import division
import numpy as np
import scipy.special as spspec

def Normal(x, y, mu, sig, corr):
    w = 1/(2*np.pi*sig[0]*sig[1]*np.sqrt(1-corr**2)) * np.exp(-1/(2*(1-corr**2)) \
    * ((x-mu[0])**2 / sig[0]**2 + (y-mu[1])**2 / sig[1]**2 - 2*corr*(x-mu[0]) \
    * (y-mu[1]) / (sig[0]*sig[1])))
    return w


def Laplace(x, y, mu, sig, corr):
    #translate sig into internal covariance structure
    Sigma = np.array([[sig[0]**2, corr*sig[0]*sig[1]], \
    [corr*sig[0]*sig[1], sig[1]**2]])
    #Force det(Gamma)==1 and Sigma== lamb*Gamma
    lamb = np.sqrt(np.linalg.det(Sigma))
    Gamma = Sigma/lamb
    
    def q(x, y, mu, Gamma):
        GammaInv = np.array([[Gamma[1,1], -Gamma[0,1]], [-Gamma[1,0], Gamma[0,0]]])
        xmu = x-mu[0]
        ymu = y-mu[1]
        qval = xmu*(GammaInv[0,0]*xmu + GammaInv[1,0]*ymu) \
        + ymu*(GammaInv[0,1]*xmu + GammaInv[1,1]*ymu)
        return qval
        
    qval = q(x, y, mu, Gamma)
    z = spspec.k0(np.sqrt(2*qval/lamb))/(np.pi*lamb)
    return z