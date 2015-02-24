# -*- coding: utf-8 -*-
"""Main solver for Epona, with supporting functions
    Written/optimized in Cython
Created on Thu May 31 20:30:27 2012

@author: Christopher Strickland
"""
#to build, run "python setup.py build_ext --inplace"
from __future__ import division
import Feedback
import numpy as np
cimport numpy as np
import time
import sys
import gc
#import odespy
import scipy.integrate
from scipy import signal
DTYPE = np.double
ctypedef np.double_t DTYPE_t

#max and min functions
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

#natural log
cdef extern from "math.h":
    double log(double)
#float epsilon
cdef double EPSILON = sys.float_info.epsilon

#treat overflow, division by zero, and invalid numbers as errors
np.seterr(all='raise',under='ignore')

cdef gfunc(np.ndarray[DTYPE_t, ndim=2] u,\
          np.ndarray[DTYPE_t, ndim=2] suit, double K, double negval):
    """Function to approximate E[Y(x)>0|Y(x_0)=0]
    
    Arguments:
        - u -- current solution, ndarray (2D)
        - suit -- ndarray, matching dimensions of u
        - K -- carrying capacity, int
    Returns:
        - g -- Approximation, ndarray"""
    cdef unsigned int rowlen, collen, row, col
    rowlen = u.shape[0]
    collen = u.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] g = np.zeros([rowlen, collen], dtype=DTYPE)
    for row in range(rowlen):
        for col in range(collen):
            if u[row,col] != 0 and suit[row,col] > negval:
                if 0 <= u[row,col] < 1:
                    #gamma correction applied in MainEqn
                    g[row,col] = -log(1-u[row,col])
                elif 1.05 >= u[row,col] >= 1:
                    #any adjustments to K here should go in function call
                    g[row,col] = K*suit[row,col]
                elif -0.01 <= u[row,col] < 0:
                    g[row,col] = 0
                else:
                    print "row "+str(row)+", col "+str(col)+", u "+str(u[row,col])
                    raise RuntimeError('u out of bounds')
    return g

cdef MainEqn(np.ndarray[DTYPE_t, ndim=2] u, double rate,\
             double gam, np.ndarray[DTYPE_t, ndim=2] suit,\
             np.ndarray[DTYPE_t, ndim=2] weightmat, int K, double negval):
    """Form of main equation
    
    Arguments:
        - rate -- Malthusian rate, float
        - gam -- Weight dependent correction value, float
        - suit -- Suitability data, ndarray
        - weightmat -- Weight matrix, ndarray
        - K -- Carrying capacity, int
        - negval -- Value under which suitability == 0, float
    Returns:
        - dudt -- Rate of change of u, ndarray"""
    cdef np.ndarray[DTYPE_t, ndim=2] dudt
    dudt = np.array(rate*gam*suit*\
        signal.fftconvolve(gfunc(u,suit,K*0.5,negval),weightmat,'same')*(1-u),\
        dtype=DTYPE)
    return dudt
    
cdef MainEqnC(np.ndarray[DTYPE_t, ndim=2] u, double rate,\
             double gam, np.ndarray[DTYPE_t, ndim=2] suit,\
             np.ndarray[DTYPE_t, ndim=2] weightmat, int K, double negval,\
             np.ndarray[DTYPE_t, ndim=2] cntrl):
    """Form of main equation with control
    
    Arguments:
        - rate -- Malthusian rate, float
        - gam -- Weight dependent correction value, float
        - suit -- Suitability data, ndarray
        - weightmat -- Weight matrix, ndarray
        - K -- Carrying capacity, int
        - negval -- Value under which suitability == 0, float
        - cntrl -- 0<= residual control <= 1 at all locations, ndarray
    Returns:
        - dudt -- Rate of change of u, ndarray"""
    cdef np.ndarray[DTYPE_t, ndim=2] dudt
    dudt = np.array(rate*gam*suit*\
        signal.fftconvolve(gfunc(u,suit,K*0.5,negval),weightmat,'same')*\
        (1-u)*(1-cntrl), dtype=DTYPE)
    return dudt
    
cdef MainEqnG(np.ndarray[DTYPE_t, ndim=2] u, double rate,\
             double gam, np.ndarray[DTYPE_t, ndim=2] suit,\
             np.ndarray[DTYPE_t, ndim=2] weightmat, int K, double negval, L,\
             G_nodeloc, G_ctrcell):
    """Form of main equation with graph
    
    Arguments:
        - rate -- Malthusian rate, float
        - gam -- Weight dependent correction value, float
        - suit -- Suitability data, ndarray
        - weightmat -- Weight matrix, ndarray
        - K -- Carrying capacity, int
        - negval -- Value under which suitability == 0, float
        - G_L -- Latent seed matrix, ndarray
        - G_nodeloc -- Node locations, list of tuples or similar
        - G_ctrcell -- center cell information
    Returns:
        - dudt -- Rate of change of u, ndarray"""
    cdef np.ndarray[DTYPE_t, ndim=2] dudt
    dudt = np.array(rate*gam*suit*\
        signal.fftconvolve(gfunc(u,suit,K*0.5,negval),weightmat,'same')*(1-u),\
        dtype=DTYPE)
    for n,pair in enumerate(G_nodeloc):
        #edge control
        for row in xrange(L.shape[1]):
            rowpos = pair[0]+row-G_ctrcell
            if 0 <= rowpos <= suit.shape[0]-1:
                for col in xrange(L.shape[2]):
                    colpos = pair[1]+col-G_ctrcell
                    if 0 <= colpos <= suit.shape[1]-1:
                        dudt[rowpos,colpos] = dudt[rowpos,colpos] +\
                        L[n,row,col]*suit[rowpos,colpos]*(1-u[rowpos,colpos])
    return dudt

cdef MainEqnGC(np.ndarray[DTYPE_t, ndim=2] u, double rate,\
             double gam, np.ndarray[DTYPE_t, ndim=2] suit,\
             np.ndarray[DTYPE_t, ndim=2] weightmat, int K, double negval, L,\
             G_nodeloc, G_ctrcell, np.ndarray[DTYPE_t, ndim=2] cntrl):
    """Form of main equation with graph and control
    
    Arguments:
        - rate -- Malthusian rate, float
        - gam -- Weight dependent correction value, float
        - suit -- Suitability data, ndarray
        - weightmat -- Weight matrix, ndarray
        - K -- Carrying capacity, int
        - negval -- Value under which suitability == 0, float
        - G_L -- Latent seed matrix, ndarray
        - G_nodeloc -- Node locations, list of tuples or similar
        - G_ctrcell -- center cell information
        - cntrl -- 0<= residual control <= 1 at all locations, ndarray
    Returns:
        - dudt -- Rate of change of u, ndarray"""
    cdef np.ndarray[DTYPE_t, ndim=2] dudt
    dudt = np.array(rate*gam*suit*\
        signal.fftconvolve(gfunc(u,suit,K*0.5,negval),weightmat,'same')*\
        (1-u)*(1-cntrl), dtype=DTYPE)
    for n,pair in enumerate(G_nodeloc):
        #edge control
        for row in xrange(L.shape[1]):
            rowpos = pair[0]+row-G_ctrcell
            if 0 <= rowpos <= suit.shape[0]-1:
                for col in xrange(L.shape[2]):
                    colpos = pair[1]+col-G_ctrcell
                    if 0 <= colpos <= suit.shape[1]-1:
                        dudt[rowpos,colpos] = dudt[rowpos,colpos] +\
                        L[n,row,col]*suit[rowpos,colpos]*(1-u[rowpos,colpos])*\
                        (1-cntrl[rowpos,colpos])
    return dudt

    
cdef CheckSubstep(np.ndarray[DTYPE_t, ndim=2] uout):
    """Check substep for values > 1, and return any correction value
    
    Arguments:
        - subu -- current solution, ndarray
        - k -- substep solution, ndarray
    Returns:
        - corval -- correction value, float.  0 if no correction needed.
        - k -- possibly fixed k value, ndarray"""
    cdef double corval
    cdef int loc
    loc = uout.argmax()
    if uout.flat[loc] > 1:
        for num,row in enumerate(uout):
            while row.max() > 1:
                uout[num,row.argmax()] = 1
        # announce clip with a -1 value
        corval = -1
    # if there is no problem, return 0
    elif np.isnan(uout.flat[loc]):
        raise RuntimeError('NaN solution detected.')
    else:
        corval = 0
    return (corval,uout)
    
cdef ArgCheck(np.ndarray[DTYPE_t, ndim=2] usub,\
              np.ndarray[DTYPE_t, ndim=2] arg):
    """Check argument of upcoming step for values > 1.
    
    Arguments:
        - subu -- current solution, ndarray
        - arg -- argument to pass to Runge Kutta substep, ndarray
    Returns:
        - corval -- correction value, float.  0 if no correction needed"""
    cdef double corval
    corval = 0
    for loc in xrange(arg.size):
        if arg.flat[loc] > 1:
            #see how big a step was taken
            stepsize = arg.flat[loc] - usub.flat[loc]
            #get smallest fraction by which to decrease h.
            if corval <= 0:
                corval = (1-usub.flat[loc])/stepsize
            else:
                corval = min(corval,(1-usub.flat[loc])/stepsize)
        elif np.isnan(arg.flat[loc]):
            raise RuntimeError('NaN solution detected.')
    return corval
    
def OptCntrl(Amat,beta,K,nm=1):
    """Calculate the optimal constant graph control.
    
    Arguments:
        - Amat -- Matrix (np.matrix)
        - beta -- beta vector (ndarray)
        - K -- scalar bound on the norm of the control (float)
        - nm -- (optional) type of norm in format of np.linalg.norm"""
    N = beta.size
    #find the kernel of Amat
    kerA = null(Amat)
    #verify that the null space is 1D
    if kerA.shape[0] != 1:
        raise ValueError('dim(ker(A)) != 1. Check that sum(nu)=1')
    #Erase negligable values
    for entry in xrange(kerA.shape[1]):
        if abs(kerA[0,entry]) < 1e-10:
            kerA[0,entry] = 0
    if np.min(kerA) < 0:
        #expect that all elements of kerA are <= or >= to zero
        kerA = -kerA
    #compute optimal control
    utild = np.zeros(beta.shape)
    #check for graph not currently relevant
    if np.linalg.norm(beta,nm) == 0:
        return utild
    for node in xrange(N):
        if kerA[0,N+node] != 0:
            utild[node] = beta[node]*(kerA[0,N+node]+kerA[0,node])/kerA[0,N+node]
    u = K*utild/np.linalg.norm(utild,nm)
    return u

def null(A, eps=1e-14):
    """Calculate the nullspace of a matrix object.
    
    Arguments:
        - A -- Matrix (np.matrix)
        - eps -- optional, largest eigenvalue to treat as zero (float)
        
    Returns:
        - null_space -- array of row vectors that span the nullspace"""
    
    u, s, v = np.linalg.svd(A)
    null_mask = (s <= eps)
    if null_mask.sum() == 0:
        null_space = np.zeros(A.shape[0])
        return np.array(null_space)
    else:
        null_space = np.compress(null_mask,v,axis=0)
        return np.array(null_space)

def GraphODEs(t,x,wa,wb,nnodes,G_nodeloc,G_ctrcell,G_days,G_Hmat,G_beta,\
        G_op,G_mu,G_Gmat,G_r,G_weightmat,G_sig,suit,usub,G_delta):
    """Definition of Graph ODEs. Assume all graph information is
        available in outer scope, and beta has been computed."""
    #S = x[:nnodes]
    #C = x[nnodes:2*nnodes]
    L = x[2*nnodes:].reshape(nnodes,wa,wb)
    
    cdef int rowpos, colpos, n, row, col, rowrange, colrange
    
    dxdt = np.concatenate((G_days*(G_Hmat.dot(x[:nnodes])-\
        G_beta*x[:nnodes]+G_op.dot(x[nnodes:nnodes*2])),\
        G_days*(G_beta*x[:nnodes]+\
        G_Gmat.dot(x[nnodes:nnodes*2])-G_mu*x[nnodes:nnodes*2])))
        
    cdef np.ndarray[DTYPE_t, ndim=3] dLdt = np.zeros((nnodes,wa,wb), dtype=DTYPE)
    conctmat = G_days*G_r*G_weightmat
    rowrange = wa
    colrange = wb
        
    #dLdt calculation
    for n,pair in enumerate(G_nodeloc):
        #edge control
        for row in xrange(rowrange):
            rowpos = pair[0]+row-G_ctrcell
            if 0 <= rowpos <= suit.shape[0]-1:
                for col in xrange(colrange):
                    colpos = pair[1]+col-G_ctrcell
                    if 0 <= colpos <= suit.shape[1]-1:
                    #combine initial spread, but cut off addition at suitability
                        dLdt[n,row,col] =\
                            conctmat[row,col]*x[nnodes+n] -\
                            G_sig*L[n,row,col]*suit[rowpos,colpos]*\
                            (1-usub[rowpos,colpos])-G_delta*L[n,row,col]
    return np.concatenate((dxdt,dLdt.flatten()))
    
def GraphODEs_Cntrl(t,x,wa,wb,nnodes,G_nodeloc,G_ctrcell,G_days,G_Hmat,G_beta,\
        G_op,G_mu,G_Gmat,G_r,G_weightmat,G_sig,suit,usub,G_delta,G_u):
    """Definition of Graph ODEs. Assume all graph information is
        available in outer scope, and beta has been computed."""
    #S = x[:nnodes]
    #C = x[nnodes:2*nnodes]
    L = x[2*nnodes:].reshape(nnodes,wa,wb)
    
    cdef int rowpos, colpos, n, row, col, rowrange, colrange
        
    dxdt = np.concatenate((G_days*(G_Hmat.dot(x[:nnodes])-\
        G_beta*x[:nnodes]+G_op.dot(x[nnodes:nnodes*2])+\
        G_u*x[nnodes:nnodes*2]),G_days*(G_beta*x[:nnodes]+\
        G_Gmat.dot(x[nnodes:nnodes*2])-(G_mu+G_u)*x[nnodes:nnodes*2])))
        
    cdef np.ndarray[DTYPE_t, ndim=3] dLdt = np.zeros((nnodes,wa,wb), dtype=DTYPE)
    conctmat = G_days*G_r*G_weightmat
    rowrange = wa
    colrange = wb
        
    #dLdt calculation
    for n,pair in enumerate(G_nodeloc):
        #edge control
        for row in xrange(rowrange):
            rowpos = pair[0]+row-G_ctrcell
            if 0<= rowpos <= suit.shape[0]-1:
                for col in xrange(colrange):
                    colpos = pair[1]+col-G_ctrcell
                    if 0 <= colpos <= suit.shape[1]-1:
                        #combine initial spread, but cut off addition at suitability
                        dLdt[n,row,col] =\
                            conctmat[row,col]*x[nnodes+n] -\
                            G_sig*L[n,row,col]*suit[rowpos,colpos]*\
                            (1-usub[rowpos,colpos])-G_delta*L[n,row,col]
    return np.concatenate((dxdt,dLdt.flatten()))
    
def Run(tspan,u0,rate,gam,suit,weightmat,K,negval,eps,cntrldict,Graph_dict=None):
    """Solve presence probabiliy equation using Runge Kutta
    
    Optimized for presence probability problem.  Variable time-step
    
    Arguments:
        - tspan -- time points to return solution at, ndarray
        - u0 -- initial presence probability, ndarray
        - rate -- Malthusian rate, float
        - gam -- Weight dependent correction value, float
        - suit -- Suitability data, ndarray
        - weightmat -- Weight matrix, ndarray
        - K -- Carrying capacity, int
        - eps -- Error bound for Runge Kutta, float
        - cntrldict -- Control parameters, dict (or None)
        - Graph_dict -- If present, run the graph coupled model
    Returns:
        - u -- Presence probability solution, ndarray (3D)
        - graph -- (If running w/ graph) Graph solution, ndarray
        - tspanout -- Time steps recorded, ndarray"""
    #Initialize solution
    cdef int a, b
    a,b = u0.shape
    u = np.zeros((tspan.size,a,b))
    u[0,:,:] = u0
    tspanout = np.zeros(tspan.shape)
    #Initialize subsolution
    usub = np.array(u0)
    #Start time at 0
    t = 0.
    #Last time step was at time 0
    tlast = 0.
    hcntr = 0
    
    #Set up graph if present
    if Graph_dict != None:
        G_days = Graph_dict['days']
        G_delta = Graph_dict['delta']
        G_gam = Graph_dict['gamma']
        G_Gmat = np.array(Graph_dict['matrix'])
        G_mu = np.array(Graph_dict['mu'])
        G_nu = np.array(Graph_dict['nu'])
        G_op = np.outer(G_nu,G_mu)
        G_Hmat = G_Gmat - np.diag(G_mu) + np.outer(G_nu,G_mu)
        G_r = Graph_dict['r']
        G_sig = Graph_dict['sigma']
        G_N = Graph_dict['N']
        G_weightmat = np.array(Graph_dict['weightmat'])
        G_ctrcell = int(G_weightmat.shape[0]/2)
        G_nodeloc = Graph_dict['nodeloc']
        G_cntrl = Graph_dict['cntrl']
        G_norm = Graph_dict['norm']
        G_K = Graph_dict['K']
        #Initialize graph solution
        nnodes = G_mu.size
        G_S = G_N*G_nu
        G_C = np.zeros(nnodes) #worth noting - we are not resetting every year
        #Initialize latent seed solution
        wa,wb = G_weightmat.shape
        G_L = np.zeros((nnodes,wa,wb))
        for num, loc in enumerate(G_nodeloc):
            rowmin = max(loc[0]-wa//2,0)
            rmindif = rowmin - (loc[0]-wa//2)
            rowmax = min(loc[0]+wa//2+1,a)
            rmaxdif = rowmax - (loc[0]+wa//2+1)
            colmin = max(loc[1]-wb//2,0)
            cmindif = colmin - (loc[1]-wb//2)
            colmax = min(loc[1]+wb//2+1,b)
            cmaxdif = colmax - (loc[1]+wb//2+1)
            if rmaxdif == 0 and cmaxdif == 0:
                G_L[num,rmindif:,cmindif:] = \
                np.array(Graph_dict['L0'])[rowmin:rowmax,colmin:colmax]
            elif rmaxdif == 0:
                G_L[num,rmindif:,cmindif:cmaxdif] = \
                np.array(Graph_dict['L0'])[rowmin:rowmax,colmin:colmax]
            elif cmaxdif == 0:
                G_L[num,rmindif:rmaxdif,cmindif:] = \
                np.array(Graph_dict['L0'])[rowmin:rowmax,colmin:colmax]
            else:
                G_L[num,rmindif:rmaxdif,cmindif:cmaxdif] = \
                np.array(Graph_dict['L0'])[rowmin:rowmax,colmin:colmax]
        #record graph states every year
        S_rec = np.zeros((tspan.size-1,nnodes))
        C_rec = np.zeros((tspan.size-1,nnodes))
        L_rec = np.zeros((tspan.size-1,nnodes,wa,wb))
        if G_cntrl != 0:
            Gu_rec = np.zeros((tspan.size-1,nnodes))
        else:
            Gu_rec = None
            
    #Set up control if present
    if cntrldict != None:
        cntrl_K = cntrldict['cntrl_K']
        cntrl_r = cntrldict['cntrl_r']
        cntrl_Lr = cntrldict['cntrl_Lr']
        import Cntrl
    
    #Starting step size
    h = 0.1
    #Old step size
    h_old = 0.
    #Last successful step size
    h_suc = 0.
    #Initialize clock
    tic = time.clock()
    #C initializations
    cdef int rowpos, colpos, row, col, n

    ### ------MAIN LOOP------ ###
    while t < tspan[-1]:
        
        #Apply control
        if cntrldict != None:
            if Graph_dict != None:
                cntrl_now, cntrlw = Cntrl.get_cntrl(\
                    usub,suit,t,tlast,cntrl_K,G_nodeloc,G_weightmat,G_C)
            else:
                cntrl_now, cntrlw = Cntrl.get_cntrl(usub,suit,t,tlast,cntrl_K)
            #if control hasn't happened yet, set to zero
            if cntrlw == None:
                cntrlw = np.zeros(suit.shape)
            if cntrl_now != None:
                usub = (1-cntrl_r*cntrl_now)*usub
                #if graph present, reduce L too
                if Graph_dict != None:
                    for n,pair in enumerate(G_nodeloc):
                        #edge control
                        for row in xrange(G_L.shape[1]):
                            rowpos = pair[0]+row-G_ctrcell
                            if 0 <= rowpos <= suit.shape[0]-1:
                                for col in xrange(G_L.shape[2]):
                                    colpos = pair[1]+col-G_ctrcell
                                    if 0 <= colpos <= suit.shape[1]-1:
                                        G_L[n,row,col] =\
                                        (1-cntrl_now[rowpos,colpos]*cntrl_Lr)*\
                                        G_L[n,row,col]
                    
        #Get correct ODEs if Graph Section is not handling it
        if Graph_dict == None or t == 0.:
            if cntrldict == None :
                odefunc = lambda x: MainEqn(x, rate, gam, suit, weightmat, K,\
                                    negval)
            else:
                odefunc = lambda x: MainEqnC(x, rate, gam, suit, weightmat, K,\
                                    negval, cntrlw)                
                                        
        #Solve underlying model for one time step of step-size h
        Feedback.FeedbackProc(0,t,h)
        h_old = h
        uout,h = RK4(usub,odefunc,h)
        Feedback.FeedbackProc(6,t,h)
        
        #GRAPH SECTION
        #this must be done after solver b/c h needs to be known
        if Graph_dict != None:
            #Calculate beta
            G_beta = np.zeros(nnodes)
            for n,pair in enumerate(G_nodeloc):
                #edge control
                for row in xrange(G_weightmat.shape[0]):
                    rowpos = pair[0]+row-G_ctrcell
                    if 0 <= rowpos <= a-1:
                        for col in xrange(G_weightmat.shape[1]):
                            colpos = pair[1]+col-G_ctrcell
                            if 0 <= colpos <= b-1:
                                #Calculate beta at each node to catch up to solution
                                G_beta[n] = G_beta[n] + G_gam*usub[rowpos,colpos]*\
                                    suit[rowpos,colpos]*G_weightmat[row,col]
                                
            #G_beta calculated. If using optimal control, calculate u if it is
            # the beginning of a year.
            if (t == 0. or (t >= tspan).sum() != (t+h >= tspan).sum())\
            and G_cntrl != 0:
                #Form A matrix.
                A = np.bmat([[G_Hmat - np.diag(G_beta), np.outer(G_nu,G_mu)],\
                             [np.diag(G_beta), G_Gmat - np.diag(G_mu)]])
                #Calculate control u
                G_u = OptCntrl(A,G_beta,G_K,G_norm)

            #Call ODE solver for h time step.
            try:
                if G_cntrl == 0:
                    #no control
                    solver = scipy.integrate.ode(GraphODEs).set_integrator(\
                        'dopri5',rtol=1e-03,atol=1e-04)
                    solver.set_f_params(wa,wb,nnodes,G_nodeloc,G_ctrcell,G_days,\
                        G_Hmat,G_beta,G_op,G_mu,G_Gmat,G_r,G_weightmat,G_sig,\
                        suit,usub,G_delta)
                else:
                    #control
                    solver = scipy.integrate.ode(GraphODEs_Cntrl).set_integrator(\
                        'dopri5',rtol=1e-03,atol=1e-04)
                    solver.set_f_params(wa,wb,nnodes,G_nodeloc,G_ctrcell,G_days,\
                        G_Hmat,G_beta,G_op,G_mu,G_Gmat,G_r,G_weightmat,G_sig,\
                        suit,usub,G_delta,G_u)
                solver.set_initial_value(\
                    np.concatenate((G_S.squeeze(),G_C.squeeze(),G_L.flatten())),t)
                solver.integrate(solver.t+h)
            except:
                import pdb, traceback
                typ, val, tb = sys.exc_info()
                traceback.print_exc()
                pdb.post_mortem(tb)
            
            #Check for negative solutions
            for n,yval in enumerate(solver.y):
                if -0.5 < yval < 0:
                    solver.y[n] = 0
            if np.sum(solver.y[:nnodes*2]) < (0.5*G_N) or\
            np.min(solver.y[:nnodes]) < 0:
                raise RuntimeError('Graph vectors seriously degraded')
            #unpack solution
            G_S = np.array(solver.y[:nnodes])
            G_C = np.array(solver.y[nnodes:nnodes*2])
            G_L = np.array(solver.y[nnodes*2:].reshape(nnodes,wa,wb))
            del solver #fix memory leak
            gc.collect()
            #renormalize solution
            G_S = G_S - (np.sum((G_S,G_C)) - G_N)*G_S/np.sum(G_S)
            #update ODE system and optimize G_L structure       
            if cntrldict == None:
                odefunc = lambda x: MainEqnG(x, rate, gam, suit, weightmat, K,\
                                      negval, G_sig*G_L, G_nodeloc, G_ctrcell)
            else:
                odefunc = lambda x: MainEqnGC(x, rate, gam, suit, weightmat, K,\
                                negval, G_sig*G_L, G_nodeloc, G_ctrcell, cntrlw)
            #record graph state every half year
            if round(t) < round(t+h):
                S_rec[int(t),:] = G_S
                C_rec[int(t),:] = G_C
                L_rec[int(t),:,:,:] = G_L
                if G_cntrl != 0:
                    Gu_rec[int(t),:] = G_u
        #END GRAPH SECTION
        
        #check for recordable step
        if (t >= tspan).sum() != (t+h >= tspan).sum():
            tdist = np.min(np.abs(t-tspan))
            thdist = np.min(np.abs(t+h-tspan))
            if tdist < thdist:
                timestep = np.argmin(np.abs(t-tspan))
                u[timestep,:,:] = np.array(usub)
                tspanout[timestep] = t
            else:
                timestep = np.argmin(np.abs(t+h-tspan))
                u[timestep,:,:] = np.array(uout)
                tspanout[timestep] = t+h
            toc = time.clock()
            eltime = toc-tic
            Feedback.ReportStep(timestep,eltime)
            
        #replace current solution with the new one
        usub = np.array(uout)
        #iterate time step
        tlast = t
        t = t+h
        #get a better stepsize
        if h_old == h:
            #This step size is not the result of internal correction
            if h != h_suc*0.945:
                #If it's a new such success, raise it a bit.  Record the old.
                h_suc = h
                h = h + h*0.05
            else:
                #stay here for a bit, raising it doesn't work
                hcntr += 1
                if hcntr == 10:
                    #we've been here for a while, raise it again
                    hcntr = 0
                    h_suc = h
                    h = h + h*0.05
        else:
            #This step size is the result of internal correction
            #Something went wrong.  Can we correct with our last success?
            if h_old != h_suc*0.945 and h_suc != 0:
                h = h_suc*0.945 #decrease a bit for safety
                #start counter for this step size, otherwise we can get stuck.
                hcntr = 0
            else:
                #the old success no longer works.
                h_suc = 0
                #leave the new h alone and see if it will go through again.
        
    ### ------ main loop finished
    Feedback.FeedbackProc(-1)
    if Graph_dict != None:
        return (u, tspanout, S_rec, C_rec, L_rec, Gu_rec)
    else:
        return (u, tspanout)
    
cdef RK4(np.ndarray[DTYPE_t, ndim=2] usub, odefunc, double h):
    """Take a step with 4th order Runge Kutta
    
    Optimized for presence probability problem.
    Arguments:
        - usub -- current solution, ndarray
        - odefunc -- RHS function handle
        - h -- current step size, float
    Returns:
        - uout -- new solution, ndarray
        - h -- step size, double"""
    ###Initialize interior solutions
    cdef unsigned int rowshp, colshp
    rowshp = usub.shape[0]
    colshp = usub.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] k1 = np.zeros([rowshp, colshp], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] k2 = np.zeros([rowshp, colshp], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] k3 = np.zeros([rowshp, colshp], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] k4 = np.zeros([rowshp, colshp], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] arg = np.zeros([rowshp, colshp], dtype=DTYPE)
    cdef double corval
    cdef np.ndarray[DTYPE_t, ndim=2] uout = np.zeros([rowshp, colshp], dtype=DTYPE)
    cdef unsigned int SFLAG, subFLAG
    
    ###Calculate interior solutions
    SFLAG = 0
    while SFLAG == 0: #test for successful solution before exit
        subFLAG = 0
        while subFLAG == 0: #finish with the exit pass, subFLAG == 1
            while True: #breaking this loop sends you back to the beginning!
                k1 = h*odefunc(usub)
                
                Feedback.FeedbackProc(1, SOLVER=4)
                
                #since the 3rd time we have usub+k3, best to overestimate now,
                # rather than run to that point.
                corval = ArgCheck(usub,usub + k1)
                if corval > 0:
                    h = h*0.95*corval #a bit smaller for safety
                    Feedback.FeedbackProc(-2, h=h)
                    break
                
                arg = usub + 0.5*k1
                k2 = h*odefunc(arg)
                Feedback.FeedbackProc(2, SOLVER=4)
                
                arg = usub + 0.5*k2
                corval = ArgCheck(usub,arg)
                if corval > 0:
                    h = h*0.95*corval #a bit smaller for safety
                    Feedback.FeedbackProc(-2,h=h)
                    break
                
                k3 = h*odefunc(arg)
                Feedback.FeedbackProc(3, SOLVER=4)
                
                arg = usub + k3
                corval = ArgCheck(usub,arg)
                if corval > 0:
                    h = h*0.95*corval #a bit smaller for safety
                    Feedback.FeedbackProc(-2,h=h)
                    break
                
                k4 = h*odefunc(arg)
                subFLAG = 1
                break
        ###Add them all up for the 4th order solution
        uout = usub + 1./6.*(k1 + 2*k2 + 2*k3 + k4)
        
        ###Check for reasonable solution
        (corval,uout) = CheckSubstep(uout)
        if corval == 0:
            SFLAG = 1 #finish
        elif corval == -1:
            Feedback.FeedbackProc(-3)
            SFLAG = 1 #finish
        else:
            h = h*0.95*corval #a bit smaller for safety
            #redo
    
    return (uout,h)