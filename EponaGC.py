# -*- coding: utf-8 -*-
"""Main Epona running file

Created on Tue Mar 06 17:22:32 2012

@author: Christopher Strickland
"""
from __future__ import division
import numpy as np
#import scipy as sp
import time
import cPickle
import DataProc
import DataIO
import EponaGSolverC
import WeightFunc2D as w2D
import Plot
import sys

def main():
    ######Prompt user to run model or plot previously generated data
    while True:
        prmtstr = ''.join(["P -- Plot previous solution\n", \
            "G -- Plot previous graph states\nQ -- Quit\n", \
            "Or press return to run the model:"])
        runmode = raw_input(prmtstr)
        runmode = runmode.strip()
        if runmode == '':
            break
        elif runmode == 'q' or runmode == 'Q':
            return
        elif runmode == "g" or runmode == "G":
            while True:
                filename = raw_input(''.join(["Please input an Epona .dat ", \
                "graph file to load or q to quit:"]))
                filename = filename.strip()
                if filename == 'q' or filename == 'Q':
                    break
                else:
                    if '_graph' not in filename:
                        filename = filename+'_graph'
                    if '.dat' not in filename:
                        filename = filename+'.dat'
                    #read data
                    try:
                        fobj = open(filename, 'r')
                        S = cPickle.load(fobj)
                        C = cPickle.load(fobj)
                        L = cPickle.load(fobj)
                        Gu = cPickle.load(fobj)
                        fobj.close()
                    except IOError as e:
                        print "I/O error: {0}".format(e.strerror)
                        continue
                    Plot.Graph(S,C,L,Gu)
                    break
        elif runmode == "p" or runmode == "P":
            while True:
                filename = raw_input(''.join(["Please input an Epona .dat ",\
                "file to load, or q to quit:"]))
                filename = filename.strip()
                if filename == 'q' or filename == 'Q':
                    break
                else:
                    if '.dat' not in filename:
                        filename = filename+'.dat'
                    #read data
                    try:
                        print('Loading data...')
                        fobj = open(filename, 'r')
                        u = cPickle.load(fobj)
                        initdata = cPickle.load(fobj)
                        runtime = cPickle.load(fobj)
                        xmesh = cPickle.load(fobj)
                        ymesh = cPickle.load(fobj)
                        bndrypoints = cPickle.load(fobj)
                        roaddata = cPickle.load(fobj)
                        presdata = cPickle.load(fobj)
                        presyears = cPickle.load(fobj)
                        nodeloc = cPickle.load(fobj)
                        fobj.close()
                    except IOError as e:
                        print "I/O error: {0}".format(e.strerror)
                        continue
                    except EOFError:
                        #Compatability with non-graph version
                        fobj.close()
                        nodeloc = None
                    #get negval
                    negval = None
                    try:
                        fobj = open('Param.ini','r')
                        for line in fobj:
                            if 'negval' in line and '=' in line:
                                (name,eq,val) = line.partition('=')
                                negval = float(val.strip())
                        fobj.close()
                    except IOError:
                        print "Param.ini not found. Using default: negval = 0.005."
                        negval = 0.005
                    except ValueError:
                        print "negval assignment corrupted in param.ini. Using default: negval = 0.005."
                    if negval is None:
                        print "negval not assigned in Param.ini.  Using default: negval = 0.005."
                        negval = 0.005
                    #plot data
                    Plot.Sol(u, initdata, runtime, xmesh, ymesh, negval,\
                        bndrypoints, roaddata, presdata, presyears,nodeloc)
            #done plotting
        else:
            print "Unrecognized input."
    #a break here corresponds to running the model.
    
    ######Load Data and Parameters
    #Prompt for graph mode
    graphmode = raw_input(\
        "Run with coupled transportation graph (Graph.ini)? [n]/y:")
    graphmode = graphmode.strip()
    if graphmode == '' or graphmode == 'n' or graphmode == 'N':
        graphmode = 0
    elif graphmode == 'y' or graphmode == 'Y':
        graphmode = 1
        import GraphIO
    else:
        print "Unrecognized input. Running without graph."
        graphmode = 0
    
    print "Processing data..."
    sys.stdout.flush()
    
    if graphmode == 1:
        # Read data and pack back into dict for passing parameters
        Graph_Param = GraphIO.ReadGraphParam()
        Graph_Param['nodes'] = DataIO.LongLatCsv(Graph_Param['graphfile'],\
            *Graph_Param['node_range'])
        Graph_Param['mu'] = GraphIO.CsvRead(Graph_Param['graphfile'],\
            *Graph_Param['out_rate_range'])
        Graph_Param['nu'] = GraphIO.CsvRead(Graph_Param['graphfile'],\
            *Graph_Param['in_prob_range'])
        #Enforce that the sum of nu = 1
        Graph_Param['nu'] = Graph_Param['nu']/np.sum(Graph_Param['nu'])
        Graph_Param['matrix'] = GraphIO.CsvRead(Graph_Param['graphfile'],\
            *Graph_Param['graph_range'])
        #Calculate correct diagonal entries
        for n,col in enumerate(Graph_Param['matrix'].T):
            Graph_Param['matrix'][n,n] =\
                -np.sum(col) + Graph_Param['matrix'][n,n]
        wargs = list(Graph_Param[k] for k in ['mean','sig','corr'])
        if Graph_Param['weight'] == 'Normal':
            G_weight = lambda x,y: w2D.Normal(x,y,*wargs)
        elif Graph_Param['weight'] == 'Laplace':
            G_weight = lambda x,y: w2D.Laplace(x,y,*wargs)
        else:
            raise ValueError('weight must be either Laplace or Normal in Param.ini')
        
    
    (Std_Param, Weight_Param) = DataIO.ReadTextParam()
    
    #Prompt for control
    if 'cntrl_K' in Std_Param:
        cntrlmode = raw_input(\
            "Run with control function (Cntrl.py)? [n]/y:")
        cntrlmode = cntrlmode.strip()
        if cntrlmode == '' or cntrlmode == 'n' or cntrlmode == 'N':
            cntrl = None
        elif cntrlmode == 'y' or cntrlmode == 'Y':
            try:
                cntrl = dict((k, Std_Param[k]) for k in \
                ('cntrl_K','cntrl_r','cntrl_Lr'))
            except KeyError as e:
                print "Parameter missing from Param.ini: {0}: ".format(e.message)
                raise
        else:
            print "Unrecognized input. Running without control function."
            cntrl = None
#    else:
#        cntrl = None
    
    (suitdata,header_suit) = DataIO.RasterAsc(Std_Param['suit_data_location'])
    initdata = DataIO.LongLatCsv(Std_Param['initial_data_location'],\
    *Std_Param['initdata_range'])
    presdata = None
    presyears = None
    if 'pres_data_location' in Std_Param:
        if Std_Param['pres_data_location'] != '':
            presdata = {}
            for n,yr in enumerate(Std_Param['pres_years']):
                datacells = Std_Param['presdata_range'][2*n:2*n+2]+\
                    Std_Param['presdata_col']
                presdata[yr] = DataIO.LongLatCsv(Std_Param['pres_data_location'],\
                    *datacells)
            presyears = Std_Param['pres_years']
            presyears.sort()
    if Std_Param['road_data_location'] != '':
        (roaddata,header_road) = DataIO.RasterAsc(Std_Param['road_data_location'])
    else:
        roaddata = None
        header_road = None
    #unpack some of the parameters for local use
    runtime = Std_Param['runtime']
    outname = Std_Param['outname']
    negval = Std_Param['negval']
        
    ######Check data for consistancy and get domain information
    (xllcorner,yllcorner,cellsize,initloc,xmesh,ymesh) = DataProc.DefineDomain(\
        initdata,suitdata,header_suit,roaddata,header_road)
        
    #Process graph nodes
    if graphmode == 1:
        Graph_Param['nodeloc'] = GraphIO.GraphProc(Graph_Param['nodes'],\
        xllcorner,yllcorner,suitdata.shape[1],suitdata.shape[0],cellsize)
        nodeloc = Graph_Param['nodes'] #for pickling
    else:
        nodeloc = None
    
    ######Process data
    #Process suitability data (by reference)
    bndrypoints = DataProc.Suit(suitdata,xllcorner,yllcorner,cellsize)   
    #Process main weight function into a spread matrix for speed
    wargs = list(Weight_Param[k] for k in ['mu','sig','corr'])
    if Weight_Param['weight'] == 'Normal':
        weight = lambda x,y: w2D.Normal(x,y,*wargs)
    elif Weight_Param['weight'] == 'Laplace':
        weight = lambda x,y: w2D.Laplace(x,y,*wargs)
    else:
        raise ValueError('weight must be either Laplace or Normal in Param.ini')
    wargs0 = list(Weight_Param[k] for k in ['mu0','sig0','corr0'])
    if Weight_Param['weight0'] == 'Normal':
        weight0 = lambda x,y: w2D.Normal(x,y,*wargs0)
    else:
        raise ValueError('internal error with initial spread weight initialization')
    arglist = list(Weight_Param[k] for k in ['weightrad','errbnd'])
    (weightmat,weighterr) = DataProc.WeightMatrix(weight, *arglist)
    
    #Process graph weight function into a spread matrix for speed
    if graphmode == 1:
        (Graph_Param['weightmat'],Graph_Param['weighterr']) =\
            DataProc.WeightMatrix(G_weight, *arglist) 
            
    #Process initial data
    u0 = DataProc.InitialPresenceData(initloc,suitdata,weight0,\
        Weight_Param['maxrad0'])
        
    #Initialize L for graph
    if graphmode == 1:
        if Graph_Param['Lnum'] != 0:
            print "Initializing latent seeds..."
            sys.stdout.flush()
            Graph_Param['L0'] = GraphIO.L0calc(u0,Graph_Param)
        else:
            Graph_Param['L0'] = np.zeros(u0.shape)
    
    ######Send information to solver
    tspan = np.arange(runtime+1)
    rate = Std_Param['rate']
    gam = Std_Param['gam']
    K = Std_Param['K']
    RKbnd = Std_Param['RKbnd']
    #Report start time
    print "Solver started on "+time.asctime()
    if graphmode == 1:
        (u,tspanout, S_rec, C_rec, L_rec, Gu_rec) = \
        EponaGSolverC.Run(tspan,u0,rate,gam,suitdata,weightmat,K,negval,RKbnd,\
        cntrl,Graph_Param)
    else:
        (u,tspanout) = \
        EponaGSolverC.Run(tspan,u0,rate,gam,suitdata,weightmat,K,negval,RKbnd,\
        cntrl)
    
    ######User feedback
    print "Solver finished on "+time.asctime()
    print "Actual time steps recorded: "+str(tspanout)
    
    ######Save data with outname
    #Pickle all data and output for local use
    fobj = open(outname+'.dat','w')
    cPickle.dump(u,fobj)
    cPickle.dump(initdata,fobj)
    cPickle.dump(runtime,fobj)
    cPickle.dump(xmesh,fobj)
    cPickle.dump(ymesh,fobj)
    cPickle.dump(bndrypoints,fobj)
    cPickle.dump(roaddata,fobj)
    cPickle.dump(presdata,fobj)
    cPickle.dump(presyears,fobj)
    cPickle.dump(nodeloc,fobj)
    fobj.close()
    print "Data has been saved for local recall in '"+outname+".dat'."
    
    #Graph states saved seperately
    if graphmode == 1:
        fobj = open(outname+'_graph.dat','w')
        cPickle.dump(S_rec,fobj)
        cPickle.dump(C_rec,fobj)
        cPickle.dump(L_rec,fobj)
        cPickle.dump(Gu_rec,fobj)
        fobj.close()
        print "Graph data has been saved for local recall in '"+outname+"_graph.dat'."
    
    fobj = open('suitdata.dat','w')
    cPickle.dump(suitdata,fobj)
    fobj.close()
    print "Suitability data saved for quick plotting recall in suitdata.dat."
    #DataIO.WriteData(u,outname)
    #print "Data has been saved to files using designation '"+outname+"'."
    
    ######Plot solutions
    Plot.Sol(u, initdata, runtime, xmesh, ymesh, negval, bndrypoints, roaddata,\
            presdata, presyears, nodeloc)


# Standard boilerplate to call the main() function to begin
# the program
if __name__ == '__main__':
    main()