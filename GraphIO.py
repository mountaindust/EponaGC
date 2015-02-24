# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Graph I/O for Epona

Created on Wed Jun 05 16:50:13 2013

@author: Christopher Strickland
"""
from __future__ import division
import numpy as np
from scipy import linalg

def null(A, eps=1e-14):
    """Calculate the nullspace of a matrix object.
    
    Arguments:
        - A -- Matrix (np.matrix)
        - eps -- optional, largest eigenvalue to treat as zero (float)
        
    Returns:
        - null_space -- array of row vectors that span the nullspace"""
    
    u, s, v = linalg.svd(A)
    null_mask = (s <= eps)
    if null_mask.sum() == 0:
        null_space = np.zeros(A.shape[0])
        return np.array(null_space)
    else:
        null_space = np.compress(null_mask,v,axis=0)
        return np.array(null_space)

def ReadGraphParam():
    """Read the Graph.ini file for graph parameters.
    
    Returns:
        - Graph_Param -- dictionary of graph parameters"""
    Params = []
    #read the data in
    try:
        fobj = open('Graph.ini','r')
        for line in fobj:
            if '=' in line and '#' not in line:
                (name,eq,val) = line.partition('=')
                val = val.strip()
                name = name.strip()
                Params.append((name,val))
            else:
                continue
        fobj.close()
    except IOError as e:
        print "I/O error: {0}".format(e.strerror)
        print "Param.ini file may be missing or corrupted."
        raise
    #process the parameters
    Param_dict = dict(Params)
    try:
        GParams = []
        GParams.append(('graphfile',Param_dict['graphfile']))
        Param_dict['node_range'] = Param_dict['node_range'].strip('()')
        GParams.append(('node_range',[int(x) for x in\
            Param_dict['node_range'].split(',')]))
        Param_dict['out_rate_range'] = Param_dict['out_rate_range'].strip('()')
        GParams.append(('out_rate_range',[int(x) for x in\
            Param_dict['out_rate_range'].split(',')]))
        Param_dict['in_prob_range'] = Param_dict['in_prob_range'].strip('()')
        GParams.append(('in_prob_range',[int(x) for x in\
            Param_dict['in_prob_range'].split(',')]))
        Param_dict['graph_range'] = Param_dict['graph_range'].strip('()')
        GParams.append(('graph_range',[int(x) for x in\
            Param_dict['graph_range'].split(',')]))
        GParams.append(('weight',Param_dict['weight'].strip()))
        Param_dict['mu'] = Param_dict['mu'].strip('()')
        GParams.append(('mean',[float(x) for x in\
            Param_dict['mu'].split(',')]))
        Param_dict['sig2'] = Param_dict['sig2'].strip('()')
        GParams.append(('sig',[float(x) for x in\
            Param_dict['sig2'].split(',')]))
        GParams.append(('corr',float(Param_dict['corr'])))
        GParams.append(('days',float(Param_dict['days'])))
        GParams.append(('delta',float(Param_dict['decay_rate'])))
        GParams.append(('gamma',float(Param_dict['gam_rate'])))
        GParams.append(('r',float(Param_dict['r_rate'])))
        GParams.append(('N',float(Param_dict['N'])))
        GParams.append(('sigma',float(Param_dict['sprout_rate'])))
        GParams.append(('Lnum',float(Param_dict['Lnum'])))
        GParams.append(('cntrl',float(Param_dict['cntrl'])))
        GParams.append(('norm',float(Param_dict['norm'])))
        GParams.append(('K',float(Param_dict['K'])))
    except ValueError as e:
        print "Parameter data error: invalid parameter\
in Graph.ini. {0}".format(e.message)
        raise
    except KeyError as e:
        print "Missing parameter in Graph.ini: {0}".format(e.message)
        raise
    #Form dictionary
    Graph_Param = dict(GParams)
    #Infer end of data cell ranges
    ndif = Graph_Param['node_range'][1]-Graph_Param['node_range'][0]
    Graph_Param['out_rate_range'].insert(1,\
        Graph_Param['out_rate_range'][0]+ndif)
    Graph_Param['out_rate_range'].append(Graph_Param['out_rate_range'][-1])
    Graph_Param['in_prob_range'].insert(1,\
        Graph_Param['in_prob_range'][0]+ndif)
    Graph_Param['in_prob_range'].append(Graph_Param['in_prob_range'][-1])
    Graph_Param['graph_range'] = Graph_Param['graph_range']+\
        [Graph_Param['graph_range'][0]+ndif,Graph_Param['graph_range'][1]+ndif]
    #switch to first row, last row, first column, last column
    Graph_Param['graph_range'][1], Graph_Param['graph_range'][2] =\
        Graph_Param['graph_range'][2], Graph_Param['graph_range'][1]
    
    #Return dictionary
    return Graph_Param
    
def CsvRead(filelocation,rowstart,rowstop,columnstart,columnstop):
    """Read a comma delimited .csv file and return the data
    
    Arguments:
        - filelocation -- string
        - rowstart -- first row to read (first row = 1), int
        - rowstop --  last row to read, int
        - columnstart -- first column to read (first col = 1), int
        - columnstop -- last column to read, int
    Returns:
        - numpy array of data"""
    import csv
    
    #Read in spreadsheet
    csvobj = csv.reader(open(filelocation,'rb'))
    npdata = np.zeros((rowstop-rowstart+1,columnstop-columnstart+1))
    for num, row in enumerate(csvobj):
        if (num >= rowstart-1) and (num < rowstop):
            npdata[num-rowstart+1,] = row[columnstart-1:columnstop]
    return npdata.squeeze()
    
def GraphProc(nodes,xllcorner,yllcorner,xlen,ylen,cellsize):
    """Convert graph nodes into (row,column) format for connecting to the
        spread model.
        
    Arguments:
        - nodes -- node locations, list of (x,y) tuples in UTM
        - xllcorner -- Easting UTM of lower left corner
        - yllcorner -- Northing UTM of lower left corner
        - cellsize -- Length of a cell edge, float
        - xlen -- Number of rows in spread model domain, int
        - ylen -- Number of columns in spread model domain, int
    Returns:
        - nodeloc -- (row,column) list of initial data locations"""
        
    nodeloc = list(0 for ii in nodes)
    for num, datapair in enumerate(nodes):
        #columns need to be done so that column 0 is highest Northing value
        nodeloc[num] = ((ylen-1)-np.floor((datapair[1] - yllcorner)/cellsize),\
            (datapair[0] - xllcorner)/cellsize)
    return nodeloc
    
def L0calc(u0,Graph_Param):
    a,b = u0.shape
    L0 = np.zeros((a,b))
    G_weightmat = Graph_Param['weightmat']
    G_nodeloc = Graph_Param['nodeloc']
    G_gam = Graph_Param['gamma']
    G_mu = Graph_Param['mu']
    G_nu = Graph_Param['nu']
    G_N = Graph_Param['N']
    G_Lnum = Graph_Param['Lnum']
    nnodes = G_mu.size
    G_ctrcell = int(G_weightmat.shape[0]/2)
    #Calculate beta
    G_beta = np.zeros(nnodes)
    for n,pair in enumerate(G_nodeloc):
        #edge control
        for row in xrange(G_weightmat.shape[0]):
            for col in xrange(G_weightmat.shape[1]):
                rowdif = row - G_ctrcell
                coldif = col - G_ctrcell
                if (pair[0]+rowdif >= 0) and (pair[0]+rowdif <= a-1)\
                and (pair[1]+coldif >= 0) and (pair[1]+coldif <= b-1):
                    #Calculate beta at each node
                    G_beta[n] = G_beta[n] +\
                        G_gam*u0[pair[0]+rowdif,pair[1]+coldif]*\
                        G_weightmat[row,col]
    #Form A matrix, find the kernel for this fixed beta
    G_Gmat = Graph_Param['matrix']
    G_Hmat = G_Gmat - np.diag(G_mu) + np.outer(G_nu,G_mu)
    A = np.bmat([[G_Hmat - np.diag(G_beta), np.outer(G_nu,G_mu)],\
                 [np.diag(G_beta), G_Gmat - np.diag(G_mu)]])
    kerA = null(A)
    if kerA.shape[0] != 1:
        raise ValueError('dim(ker(A)) != 1. Check that sum(nu)=1')
    #Erase negligable values
    for node in xrange(kerA.shape[1]):
        if abs(kerA[0,node]) < 1e-10:
            kerA[0,node] = 0
    if np.min(kerA) < 0:
        #expect that all elements of kerA are <= or >= to zero
        kerA = -kerA
    #Calculate the equilibrum solution
    xstar = G_N*kerA/linalg.norm(kerA,1)
    xstar = xstar.flatten()
    
    #Spread some seeds to L0 based on xstar
    for n,pair in enumerate(G_nodeloc):
        #edge control
        for row in xrange(G_weightmat.shape[0]):
            for col in xrange(G_weightmat.shape[1]):
                rowdif = row - G_ctrcell
                coldif = col - G_ctrcell
                if (pair[0]+rowdif >= 0) and (pair[0]+rowdif <= a-1)\
                and (pair[1]+coldif >= 0) and (pair[1]+coldif <= b-1):
                    #combine initial spread, but cut off addition at suitability
                    L0[pair[0]+rowdif,pair[1]+coldif] =\
                        G_Lnum*G_weightmat[row,col]*xstar[nnodes+n]
    return L0