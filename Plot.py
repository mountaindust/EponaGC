# -*- coding: utf-8 -*-
"""Plotting for Epona
Created on Tue Jun 19 11:32:41 2012

@author: Christopher Strickland
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcl
import matplotlib.cm as cm
import matplotlib as mpl

def Sol(u, initdata, runtime, xmesh, ymesh, negval, bndrypoints,\
        roaddata=None,presdata=None,presyears=None,nodeloc=None):
    """Prompt user for time step and plot Epona solution
    
    Arguments:
        -u
        -initdata
        -runtime
        -xmesh
        -ymesh
        -negval
        -bndrypoints
        -roadloc
        -presdata=None
        -presyear=None
        -suitability=None
    """
    #Norm definition
    mynorm = mcl.Normalize(0,1,clip=True)
    #Get bndry and road locations
    xbndry = []
    ybndry = []
    xinit = []
    yinit = []
    for row in bndrypoints:
        xbndry.append(row[0])
        ybndry.append(row[1])
    for row in initdata:
        xinit.append(row[0])
        yinit.append(row[1])
    if roaddata is not None:
        try:
            roadmask = np.ma.array(roaddata[::-1,:],\
                mask = roaddata[::-1,:] != 1, fill_value = 0)
            roadnorm = mcl.Normalize(0,1)
        except:
            print 'Error in road raster plot. Layer omitted.'
            roadmask = None
    
    #plot selected time steps
    suit_mode = -1
    while True:
        tt = raw_input("Input time step to plot, 0-"+str(runtime)+\
            ", s to use/setup suitability data, or q to quit:")
        tt = tt.strip()
        if tt == "q" or tt == "Q":
            break
        elif tt == "s" or tt == "S":
            
            if suit_mode == -1:
                #no suitability loaded
                import cPickle
                filename = raw_input("Input suitability .dat file to load, "+\
                    "or press enter to load suitdata.dat:")
                filename = filename.strip()
                if filename == '':
                    try:
                        fobj = open('suitdata.dat','r')
                        suit = cPickle.load(fobj)
                        fobj.close()
                        #flip for plotting
                        suit = suit[::-1,:]
                        suit_mode = 0
                    except IOError as e:
                        print "I/O error: {0}".format(e.strerror)
                        continue
                else:
                    if '.dat' not in filename:
                        filename = filename+'.dat'
                    #read data
                    try:
                        fobj = open(filename, 'r')
                        suit = cPickle.load(fobj)
                        fobj.close()
                        #flip for plotting
                        suit = suit[::-1,:]
                        suit_mode = 0
                    except IOError as e:
                        print "I/O error: {0}".format(e.strerror)
                        continue
                    
            #once loaded:        
            ss = raw_input("0 - do not plot suitability\n"+\
                           "1 - plot suitabilty separately\n"+\
                           "2 - scale presence solution by local suitability\n"+\
                           "3 - do both\n:")
            try:
                ss = int(ss)
            except ValueError:
                print "Unrecognized input. No change made."
                continue
            if ss > 3 or ss < 0:
                print "Unrecognized input. No change made."
                continue
            suit_mode = ss
            continue
        else:
            try:
                tt = int(tt)
            except ValueError:
                print "Unrecognized input."
                continue
        if tt < 0 or tt > runtime:
            print "Input is outside of the data range."
            continue
        if suit_mode < 2:
            #flip matrix updown-pcolor plots upsidedown
            plotvalues = np.ma.array(u[tt,::-1,:],\
                mask = u[tt,::-1,:] < negval, fill_value = 0)
        else:
            plotvalues = np.ma.array(u[tt,::-1,:]*suit,\
                mask = u[tt,::-1,:] < negval, fill_value = 0)
                
        #prompt for color or B/W
        c_switch = raw_input(''.join(["Hit enter for color ", \
                "or input any character for B/W:"]))
        c_switch = c_switch.strip()
        if c_switch == '':
            #color
            #cmap = cm.gist_earth_r
            #take a portion of the gist_earth_r colormap
            cmap = cm.colors.LinearSegmentedColormap.from_list(\
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cm.gist_earth_r,\
                a=0.15, b=0.85), cm.gist_earth_r(np.linspace(0.15,0.85,100)))
        else:
            #grayscale. want to leave out black, so need a new colormap
            #get a grayscale colormap with 10000 steps
            rgba_array = cm.binary(np.linspace(0,1,num=10000,endpoint=True))
            #get 10%-70% of it. This is our new colormap.
            extract_rgba_array_255 = rgba_array[1000:7000,0:3]
            cmap = cm.colors.ListedColormap(extract_rgba_array_255)
        
        print('Rendering...')
        cmap.set_bad('w',1.)
        mpl.rc('xtick', labelsize=18)
        mpl.rc('ytick', labelsize=18)
        plt.pcolor(xmesh,ymesh,plotvalues,norm=mynorm,cmap=cmap)
        plt.autoscale(tight=True)
        plt.colorbar(ticks=np.linspace(0,1,11))
        plt.xlabel('UTM Easting',fontsize=24)
        plt.ylabel('UTM Northing',fontsize=24)
        if suit_mode < 2:
            plt.title('Species Range Probability at Year '+str(tt),\
            fontsize=36)
        else:
            plt.title('Species Presence Probability at Year '+str(tt),\
            fontsize=36)
        plt.hold(True)
        #plot boundary
        plt.scatter(xbndry,ybndry,s=10,c='k',marker='D')
        #plot roads
        if roadmask is not None:
            try:
                cmaproad = cm.binary
                cmaproad.set_bad(alpha=0) #transparent
                plt.pcolor(xmesh,ymesh,roadmask,norm=roadnorm,cmap=cmaproad)
            except:
                print 'Error in road raster plot. Layer omitted.'
        
        #show original presence data
        if c_switch == '':
            #color
            plt.scatter(xinit,yinit,s=100,c=(0.85,0,0),marker='x',linewidths=3)
        else:
            #grayscale
            plt.scatter(xinit,yinit,s=100,c='k',marker='x',linewidths=3)
            plt.scatter(xinit,yinit,s=80,c='w',marker='x')
        
        #plot presence data for most recent year
        if presdata is not None:
            #find the last year <= than the current year in sorted list
            year = -1
            for yr in presyears:
                if yr < tt:
                    year = yr
                elif yr == tt:
                    year = yr
                    break
                else:
                    break
            if year > -1:
                xpres = []
                ypres = []
                for row in presdata[year]:
                    xpres.append(row[0])
                    ypres.append(row[1])
                if c_switch == '':
                    #color
                    plt.scatter(xpres,ypres,s=120,c='m',marker='+',linewidths=3)
                else:
                    #grayscale
                    plt.scatter(xpres,ypres,s=120,c='k',marker='+',linewidths=3)
                    plt.scatter(xpres,ypres,s=96,c='w',marker='+')
                
        #plot graph node locations
        if nodeloc is not None:
            xloc = []
            yloc = []
            for row in nodeloc:
                xloc.append(row[0])
                yloc.append(row[1])
            if c_switch == '':
                #color
                plt.scatter(xloc,yloc,s=120,c=(0,.75,0),marker='H')
            else:
                #grayscale
                plt.scatter(xloc,yloc,s=120,c='w',marker='H')
        
        plt.hold(False)      
        
        if suit_mode == 1 or suit_mode == 3:
            plt.figure()
            cmap.set_bad('w',1.)
            mpl.rc('xtick', labelsize=18)
            mpl.rc('ytick', labelsize=18)
            plt.pcolor(xmesh,ymesh,suit,norm=mynorm,cmap=cmap)
            plt.autoscale(tight=True)
            plt.colorbar(ticks=np.linspace(0,1,11))
            plt.xlabel('UTM Easting',fontsize=24)
            plt.ylabel('UTM Northing',fontsize=24)
            plt.title('Species Suitability Map',fontsize=36)
            plt.hold(True)
            #plot boundary
            plt.scatter(xbndry,ybndry,s=10,c='k',marker='D')
        plt.show()
        
def Graph(S,C,L,Gu):
    runtime, nnodes = S.shape
    
    #Graph vectors
    plt.figure()
    plt.subplot(211)
    plt.hold(True)
    for tt in range(runtime):
        plt.scatter(range(nnodes),S[tt,:],c=str(0.75-0.75*tt/runtime))
    plt.xlabel('Node number')
    plt.ylabel('Population')
    plt.title('Susceptible graph vectors')
    plt.hold(False)
    plt.subplot(212)
    plt.hold(True)
    for tt in range(runtime):
        plt.scatter(range(nnodes),C[tt,:],c=str(0.75-0.75*tt/runtime))
    plt.xlabel('Node number')
    plt.ylabel('Population')
    plt.title('Carrier graph vectors')
    plt.hold(False)
    
    #Control
    if Gu is not None:
        plt.figure()
        plt.hold(True)
        for tt in range(runtime):
            plt.scatter(range(nnodes),Gu[tt,:],c=str(0.75-0.75*tt/runtime))
        plt.xlabel('Node number')
        plt.ylabel('Control rate')
        plt.title('Graph control over time')
        plt.hold(False)
        
    #Latent seeds
    #Plot sum of latent seeds in each node location
    Lsum = np.sum(L,(2,3))
    plt.figure()
    plt.hold(True)
    for tt in range(runtime):
        plt.scatter(range(nnodes),Lsum[tt,:],c=str(0.75-0.75*tt/runtime))
    plt.xlabel('Node number')
    plt.ylabel('Total latent seeds')
    plt.title('Latent seed numbers')
    plt.hold(False)
    plt.show()