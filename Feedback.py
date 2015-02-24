# -*- coding: utf-8 -*-
"""Solver feedback
Feedback functions for EponaSolver
Created on Thu Jun 07 14:12:45 2012

@author: Christopher Strickland
"""
from __future__ import division
import sys
import time
import numpy as np

def FeedbackProc(funcdone, t=None,h=None,SOLVER=5):
    """Provide user feedback on solution progress
    
    Arguments:
        - funcdone -- step just completed, (0-6,-1 done,-2 h too large,-3 clip)
        - t -- optional, current time step, float
        - h -- optional, current step size, float
        - SOLVER -- optional, solver order for inner timesteps"""
    if funcdone == 0:
        print "Runge Kutta starting at t={0} with step size={1}.".format(t,h)
    elif funcdone == 6:
        print "done."
    elif funcdone == -1:
        print "All done!"
    elif funcdone == -2:
        print "h too large. Restarting with step size={0}...".format(h)
    elif funcdone == -3:
        sys.stdout.write("clip-")
    else:
        if SOLVER == 5:
            sys.stdout.write("{}/6...".format(funcdone))
        elif SOLVER == 4:
            sys.stdout.write("{}/4...".format(funcdone))
    sys.stdout.flush()
    
def ReportFailure(failnum, t, h):
    """Report RKCK step failure, step to be redone
    
    Arguments:
        - failnum -- failure number, int
        - t -- failed time step, float
        - h -- failed step size, float"""
    print "Runge Kutta step failure #{0} reported at t={1} with step size {2}. Repeating step.\n".format(failnum,t,h)
    
def ReportStep(number, eltime):
    """Report recorded time step
    
    Arguments:
        - number -- time step recorded, float
        - eltime -- elapsed time, float"""
    print 'Time step '+str(number)+' recorded at '+time.asctime()
    eldays = int(np.floor(eltime/60./60./24.))
    elhrs = int(np.floor((eltime - eldays*60.*60.*24.)/60./60.))
    elmin = int(np.floor((eltime - eldays*60.*60.*24. - elhrs*60.*60.)/60.))
    elsec = eltime - eldays*60.*60.*24. - elhrs*60.*60. - elmin*60.
    print 'Time elapsed: '+str(eldays)+':'+str(elhrs)+':'+str(elmin)+':'+str(elsec)