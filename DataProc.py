# -*- coding: utf-8 -*-
"""Data Processing for Epona

Created on Wed May 30 16:35:28 2012

@author: Christopher Strickland
"""
from __future__ import division
import numpy as np

def IntWeight(weight,xmin,xmax,ymin,ymax):
    """Integrate a given 2D weight function
    
    Arguments:
        - weight -- inline function handle to be integrated
        - xmin -- float
        - xmax -- float
        - ymin -- float
        - ymax -- float
    Returns:
        - val -- value of integration within rectangle"""
    gridsize = 18 #controls degree of accuracy. MUST BE AN EVEN NUMBER
                  #TO DEAL WITH INF AT ORIGIN (E.G. LAPLACE DIST.).
                  #18 should provide nearly three digits of accuracy with
                  #Laplace when sig=[1,1]
                  
    intarea = np.zeros((gridsize,gridsize))
    xstepsize = abs(ymax - ymin)/gridsize
    ystepsize = abs(xmax - xmin)/gridsize
    cellsize = xstepsize*ystepsize
    for stepx in range(0,gridsize):
        for stepy in range(0,gridsize):
            #locate yourself in the middle of a subcell
            xx = xmin + xstepsize/2. + xstepsize*stepx
            yy = ymin + ystepsize/2. + ystepsize*stepy
            #multiply the function times the area of the square
            #INCLUDE ERROR HANDLE TO DEAL WITH inf
            intarea[stepx,stepy] = weight(xx,yy)*cellsize
    #sum the total
    val = intarea.sum()
    return val

def DefineDomain(initdata,suitdata,header_suit,roaddata,header_road):
    """Get domain information, form xmesh/ymesh, and check for inconsistancies
    
    Arguments:
        - initdata -- initial data, list of (x,y) tuples
        - suitdata -- suitability data, ndarray
        - header_suit -- suitability header, dictionary
        - roaddata -- road raster data, ndarray
        - header_road -- road raster header, dictionary
    Returns:
        - xllcorner -- Easting UTM of lower left corner
        - yllcorner -- Northing UTM of lower left corner
        - cellsize -- Length of a cell edge, float
        - initloc -- (row,column) list of initial data locations
        - xmesh -- ndarray of UTM Easting coordinates
        - ymesh -- ndarray of UTM Northing coordinates"""
    
    if header_suit != header_road:
        raise Warning('Headers for suitability data and road data do not match.')
    if roaddata is not None:
        if suitdata.shape != roaddata.shape:
            raise Warning('Raster shape for suitability data and road data'+\
            ' do not match. Road layer may not plot correctly.')
        
    #Get domain information if header was present
    if header_suit:
        xllcorner = header_suit['xllcorner']
        yllcorner = header_suit['yllcorner']
        cellsize = float(header_suit['cellsize'])
    else:
        raise NameError('Undefined domain: header_suit is empty')
        
    (ylen,xlen) = suitdata.shape
    
    #Convert initdata into (row,column) format for processing into u0
    initloc = list(0 for ii in initdata)
    for num, datapair in enumerate(initdata):
        #columns need to be done so that column 0 is highest Northing value
        initloc[num] = ((ylen-1) - np.floor((datapair[1] - yllcorner)/cellsize),\
            np.floor((datapair[0] - xllcorner)/cellsize))
    
    #Create domain mesh
    #xmesh and ymesh are lines on bndry of cells.
    #Make sure to get the right most and top most coordinates!
    xmesh = np.arange(xllcorner,xllcorner + xlen*cellsize+2,cellsize)
    ymesh = np.arange(yllcorner,yllcorner + ylen*cellsize+2,cellsize)
    
    #complain if data is outside suitdata range (?)
    
    return (xllcorner,yllcorner,cellsize,initloc,xmesh,ymesh)

def Suit(suitdata,xllcorner,yllcorner,cellsize):
    """Process suitability data for boundaries and excessively low suitability
    
    Arguments:
        - suitdata -- suitability data, ndarray
        - xllcorner -- float
        - yllcorner -- float
        - cellsize -- float
    Returns 
        - suitdata -- reformatted, ndarray
        - bndrypoints -- (x,y) UTM list"""
    #initialize
    bounds = np.zeros(suitdata.shape)
    #suitdata is passed by reference. Process.
    for row in range(0,suitdata.shape[0]):
        for col in range(0,suitdata.shape[1]):
            if 1 >= suitdata[row,col] >= 0.001:
                bounds[row,col] = 1
            elif 0.001 > suitdata[row,col] >= 0:
                suitdata[row,col] = 0
                bounds[row,col] = 1
            else:
                suitdata[row,col] = 0
                bounds[row,col] = 0
    #sort through boundary raster for boundary points
    bndrypoints = []
    for row in range(1,bounds.shape[0]-1):
        for col in range(1,bounds.shape[1]-1):
            if (bounds[row-1,col] == 0 or \
            bounds[row+1,col] == 0 or bounds[row,col-1] == 0 or \
            bounds[row,col+1] == 0 ) and bounds[row,col] == 1:
                xshift = (col+0.5)*cellsize #center in cell
                yshift = (bounds.shape[0]-row-0.5)*cellsize
                bndrypoints.append((xllcorner+xshift,yllcorner+yshift))
    return bndrypoints

def RoadLoc(roaddata,xllcorner,yllcorner,cellsize):
    """Process road data for plotting. Depreciated.
    
    Assumes that road data area matches perfectly with suitdata area
    
    Arguments:
        - roaddata -- raster of road data, ndarray
        - xllcorner -- float
        - yllcorner -- float
        - cellsize -- float    
    Returns: 
        - roadloc -- (x,y) list for road locations in UTM"""
    roadloc = []
    for row in range(0,roaddata.shape[0]):
        for col in range(0,roaddata.shape[1]):
            if roaddata[row,col] == 1:
                xshift = (col+0.5)*cellsize #center in cell
                yshift = (roaddata.shape[0]-row-0.5)*cellsize
                roadloc.append((xllcorner+xshift,yllcorner+yshift))
    return roadloc
    
def WeightMatrix(weight,weightrad,errbnd):
    """Process the main weight function into a spread matrix for speed
    
    Interpret errbnd == 0 as a demand for the matrix to be constructed with
    weightrad, exactly.  In such a case, report the error.
    
    Arguments:
        - weight -- inline function handle
        - weightrad -- starting radius, int
        - errbnd -- maximum area outside matrix, float
    Returns:
        - weightmat -- weight matrix, ndarray
        - weighterr -- area left outside of matrix, redistributed evenly, float"""
    if errbnd == 0:
        sprdlen = weightrad*2 + 1
        weightmat = np.zeros((sprdlen,sprdlen))
        ctrcell = weightrad       
        for col in range(sprdlen):
            for row in range(sprdlen):
                if np.sqrt((row-ctrcell)**2 + (col-ctrcell)**2) <= weightrad:
                    #integrate more accurately each square in the matrix
                    weightmat[row,col] = IntWeight(weight,(col-ctrcell)-0.5,\
                        (col-ctrcell)+0.5,(row-ctrcell)-0.5,(row-ctrcell)+0.5)
        #report error
        weighterr = np.abs(1-weightmat.sum())
    else:
        #get an educated guess as to the radius needed
        strmesh = range(-30,31)
        centr = 30 #index of center cell in strmesh
        testx = np.zeros(np.size(strmesh))
        testy = np.zeros(np.size(strmesh))
        #integrate the cells along the x and y axis
        for tt in range(np.size(strmesh)):
            testx[tt] = IntWeight(\
                weight,strmesh[tt]-0.5,strmesh[tt]+0.5,-0.5,0.5)
            testy[tt] = IntWeight(\
                weight,-0.5,0.5,strmesh[tt]-0.5,strmesh[tt]+0.5)
        #find the maximum integration values between the two for each radius
        test = np.max(np.array([testx,testy]),0)
        #estimate upper bound for area at each radius
        for tt in range(np.size(test)):
            test[tt] = test[tt]*2*np.pi*(np.abs(strmesh[tt]+0.5))
        #enforce weightrad as the minimum radius to integrate over
        testval = np.max(np.array([test[centr-weightrad],test[centr+weightrad]]))
        #search from weightrad, increasing the radius, until a good 
        #start radius is found
        while (testval > errbnd) and (weightrad+centr < np.size(strmesh)):
            weightrad = weightrad + 1
            testval = np.max(np.array([test[centr-weightrad],test[centr+weightrad]]))
        #loop breaks when good radius is found, or radius > 30
        #a little cleanup now...
        del testval, test, testx, testy, centr, strmesh
        #integrate at the starting radius, increasing the radius until below
        #specified error bound
        while True:
            sprdlen = weightrad*2 + 1
            weightmat = np.zeros((sprdlen,sprdlen))
            ctrcell = weightrad
            for row in range(sprdlen):
                for col in range(sprdlen):
                    if np.sqrt((row-ctrcell)**2 + (col-ctrcell)**2) <= weightrad:
                        #integrate more accurately each square in the matrix
                        weightmat[row,col] = IntWeight(weight,(col-ctrcell)-0.5,\
                            (col-ctrcell)+0.5,(row-ctrcell)-0.5,(row-ctrcell)+0.5)
            #find error
            weighterr = np.abs(1-weightmat.sum())
            if weighterr <= errbnd:
                break
            else:
                weightrad = weightrad + 1
    #distribute error evenly across the weight matrix
    (a,b) = weightmat.shape
    weightmat = weightmat + weighterr/(a*b)
    
    return (weightmat,weighterr)
    
def InitialPresenceData(initloc,suitdata,weight0,maxrad0):
    """Given initial data as species presence locations, create initial conditions for presence probability
            
    Arguments:
        - initloc -- initial data, (row,column) list
        - suitdata -- suitability data, ndarray
        - weight0 -- inline function handle
        - maxrad0 -- furthest to spread initial conditions, int
    Returns:
        - u0 -- initial conditions for species presence, ndarray"""
    rowlen = suitdata.shape[0]
    collen = suitdata.shape[1]
    u0 = np.zeros((rowlen,collen))
    sprdlen0 = maxrad0*2 + 1
    sprdmat0 = np.zeros((sprdlen0,sprdlen0))
    ctrcell0 = maxrad0
    #calculate spread matrix
    for row in range(sprdlen0):
        for col in range(sprdlen0):
            if np.sqrt((row-ctrcell0)**2 + (col-ctrcell0)**2) <= maxrad0:
                sprdmat0[row,col] = IntWeight(weight0,(col-ctrcell0)-0.5,\
                    (col-ctrcell0)+0.5,(row-ctrcell0)-0.5,(row-ctrcell0)+0.5)
    #scale the matrix so that the origin == 1
    sprdmat0 = sprdmat0/sprdmat0[ctrcell0,ctrcell0]
    #spread initial conditions
    for pair in initloc:
        for row in range(sprdlen0):
            for col in range(sprdlen0):
                rowdif = row - ctrcell0
                coldif = col - ctrcell0
                if (pair[0]+rowdif >= 0) and (pair[0]+rowdif <= rowlen-1) and\
                (pair[1]+coldif >= 0) and (pair[1]+coldif <= collen-1):
                    #combine initial spread, but cut off addition at suitability
                    u0[pair[0]+rowdif,pair[1]+coldif] = np.min(np.array(\
                    [u0[pair[0]+rowdif,pair[1]+coldif] +\
                    suitdata[pair[0]+rowdif,pair[1]+coldif] * sprdmat0[row,col],\
                    suitdata[pair[0]+rowdif,pair[1]+coldif]]))
    #initial condition at presence location = 1, regardless of suitability
    for pair in initloc:
        u0[pair[0],pair[1]] = 1
    
    return u0
        
def InitialPopData(initloc,suitdata,K,weight0,maxrad0):
    """Given initial data as species population, create initial conditions
        for presence probability
        
    Arguments:
        - initloc -- initial data, (row,column) list
        - suitdata -- suitability data, ndarray
        - K -- carrying capacity, int
        - weight0 -- inline function handle
        - maxrad0 -- furthest to spread initial conditions, int
    Returns:
        - u0 --  initial conditions for species presence, ndarray"""
    pass
