# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Landscape based control module.
Edit this module to create custom control regimes which can then be loaded
into Epona.

Epona will call get_cntrl, which should be edited to call the desired function
in the module.

All control functions are expected to return relative amounts of control to be
applied.  0  is no control.  All numbers will be normalized.

Created on Thu Jul 18 13:08:28 2013

@author: Christopher Strickland
"""
from __future__ import division
import numpy as np
from math import fmod
from math import pow
#import scipy as sp

def node_cntrl(sol,suit,time,tlast,cntrl_K,node_loc,Gweightmat):
    """This function provides a control in suitable (>#) areas
    around node locations"""
    applytime = 0.1 #mod 1
    
    #control will take place every year at t=0.1
    if fmod(tlast,1) < applytime and fmod(time,1) >= applytime:
        cntrl_now = np.zeros(suit.shape)
        cntrlw = np.zeros(suit.shape)
        max_dist = Gweightmat.shape[0]//2 #max distance away from node to consider
        for node in node_loc:
            #preference locations within node range first
            for ii in xrange(-max_dist,max_dist):
                row = node[0]+ii
                ## edge control - always include this code!
                if 0 <= row < suit.shape[0]:
                    for jj in xrange(-max_dist,max_dist):
                        col = node[1]+jj
                        #edge control condition
                        if 0 <= col < suit.shape[1]:
                ## end edge control
                            #### control condition ####
                            if suit[row,col]*sol[row,col] > 0.1:
                                cntrlw[row,col] = pow(sol[row,col],0.5)
        if cntrlw.sum() < cntrl_K:
            #amp up the parts that are less than 1
            extraK = cntrl_K - cntrlw.sum()
            tempmat = np.zeros(suit.shape)
            for ii in xrange(suit.shape[0]):
                for jj in xrange(suit.shape[1]):
                    if 0 < cntrlw[ii,jj] < 1:
                        tempmat[ii,jj] = cntrlw[ii,jj]
            if (1 - tempmat).sum() <= extraK:
                for ii in xrange(suit.shape[0]):
                    for jj in xrange(suit.shape[1]):
                        if 0 < cntrlw[ii,jj] < 1:
                            cntrlw[ii,jj] = 1
            else:
                tempmat = 1 - tempmat
                for ii in xrange(suit.shape[0]):
                    for jj in xrange(suit.shape[1]):
                        if 0 < cntrlw[ii,jj] < 1:
                            cntrlw[ii,jj] = cntrlw[ii,jj] +\
                                tempmat[ii,jj]*extraK/tempmat.sum()
        elif cntrlw.sum() > cntrl_K:
            #reduce globally to meet constraint (other methods too costly)
            cntrlw = cntrlw*cntrl_K/cntrlw.sum()
        cntrl_now = np.array(cntrlw)
        node_cntrl.cntrlw = np.array(cntrlw)
    else:
        cntrl_now = None
        # Decay of lingering control
        if node_cntrl.cntrlw != None:
            # No decay for first while
            if 0 <= fmod(time,1)-applytime < 0.2:
                cntrlw = node_cntrl.cntrlw
            elif 0.2 <= fmod(time,1)-applytime:
                t = fmod(time,1)-applytime-0.2
                tend = 1-applytime-0.2
                cntrlw = node_cntrl.cntrlw - pow((t/tend),2)*node_cntrl.cntrlw
            else:
                cntrlw = np.zeros(suit.shape)
        else:
            cntrlw = node_cntrl.cntrlw
    
    return cntrl_now, cntrlw
# Initialization of persistant variable (used to prevent costly recalculation)
node_cntrl.cntrlw = None

def get_cntrl(sol,suit,time,tlast,cntrl_K,node_loc=None,Gweightmat=None,GC=None):
    """Method that provides an interface between this module and the solver.
    Hopefully prevents needing to recompile the solver again and again.
    
    Arguments:
        - sol -- Current solution, ndarray (2D)
        - suit -- suitability, ndarray (2D)
        - time -- t in solver, float
        - tlast -- t of previous time step in solver, float
        - cntrl_K -- max 1-norm amount of control, float
        - node_loc -- (optional) node locations in solution array, list of tuples
        - Gweightmat -- (optional) weight matrix connecting graph nodes to
                        physical location
        - GC -- (optional) current carrier vector solution, ndarray (1D)
        
    Returns:
        - cntrl_now -- amount of immediate kill, between 0 and 1, ndarray (2D) or None.
        - cntrlw -- residual control, between 0 and 1, ndarray (2D)
    """
    
    ### Overwrite with control function you wish you use ###
    ### The function should return an ndarray the same size as sol or suit
    ### With relative amounts of control to be applied in each location
    cntrl_now, cntrlw = node_cntrl(sol,suit,time,tlast,cntrl_K,node_loc,Gweightmat)
    ###
    
    return cntrl_now, cntrlw