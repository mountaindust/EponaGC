# -*- coding: utf-8 -*-
"""Data I/O for Epona

Created on Thu Feb 16 13:47:37 2012

@author: Christopher Strickland
"""

from __future__ import division
import numpy as np
try:
    from xlrd import open_workbook
    xlrd_FLAG = True
except ImportError:
    xlrd_FLAG = False

def ReadTextParam():
    """Read the Param.ini file for program parameters.
    
    Returns:
        - Std_Param -- dictionary of standard parameters
        - Weight_Param -- dictionary of weight parameters"""
    Params = []
    #read the data in
    try:
        fobj = open('Param.ini','r')
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
        #Standarad parameters
        SParams = []
        SParams.append(('runtime', int(Param_dict['runtime'])))
        SParams.append(('rate', float(Param_dict['rate'])))
        #check to see if K is an int
        K = int(Param_dict['K'])
        #save K as a float
        SParams.append(('K', float(K)))
        SParams.append(('gam', float(Param_dict['correction'])))
        SParams.append(('negval', float(Param_dict['negval'])))
        Param_dict['initdata_range'] = Param_dict['initdata_range'].strip('()')
        SParams.append(('initdata_range', [int(x) for x in\
            Param_dict['initdata_range'].split(',')]))
        if Param_dict['out_prefix'] == '':
            raise ValueError('out_prefix must take a value in Param.ini')
        else:
            SParams.append(('outname', Param_dict['out_prefix']))
        if Param_dict['suitability_data'] == '':
            raise ValueError('suitability_data must be designated in Param.ini')
        else:
            SParams.append(('suit_data_location', Param_dict['suitability_data']))
        if Param_dict['initial_data'] == '':
            raise ValueError('initial_data must be designated in Param.ini')
        else:
            SParams.append(('initial_data_location', Param_dict['initial_data']))
        SParams.append(('road_data_location', Param_dict['road_data']))
        if 'pres_data' in Param_dict:
            SParams.append(('pres_data_location', Param_dict['pres_data']))
            Param_dict['presdata_col'] = Param_dict['presdata_col'].strip('()')
            SParams.append(('presdata_col',[int(x) for x in\
                Param_dict['presdata_col'].split(',')]))
            Param_dict['presdata_range'] = Param_dict['presdata_range'].strip('()')
            SParams.append(('presdata_range',[int(x) for x in\
                Param_dict['presdata_range'].split(',')]))
            Param_dict['pres_years'] = Param_dict['pres_years'].strip('()')
            SParams.append(('pres_years', [int(x) for x in\
                Param_dict['pres_years'].split(',')]))
        #Weight parameters
        WParams = []
        Param_dict['weight'] = Param_dict['weight'].capitalize()
        WParams.append(('weight',Param_dict['weight']))
        Param_dict['mu'] = Param_dict['mu'].strip('()')
        WParams.append(('mu',[float(x) for x in Param_dict['mu'].split(',')]))
        Param_dict['sig2'] = Param_dict['sig2'].strip('()')
        WParams.append(('sig',\
            [np.sqrt(float(x)) for x in Param_dict['sig2'].split(',')]))
        WParams.append(('corr',float(Param_dict['corr'])))
        Param_dict['mu_0'] = Param_dict['mu_0'].strip('()')
        WParams.append(('mu0',[float(x) for x in Param_dict['mu_0'].split(',')]))
        Param_dict['sig2_0'] = Param_dict['sig2_0'].strip('()')
        WParams.append(('sig0',\
            [np.sqrt(float(x)) for x in Param_dict['sig2_0'].split(',')]))
        WParams.append(('corr0',float(Param_dict['corr_0'])))
    except ValueError as e:
        print "Parameter data error: invalid parameter\
in Param.ini. {0}".format(e.message)
        raise
    except KeyError as e:
        print "Missing parameter in Param.ini: {0}".format(e.message)
        raise
    #Control parameters
    try:
        SParams.append(('cntrl_K',float(Param_dict['cntrl_K'])))
        SParams.append(('cntrl_r',float(Param_dict['cntrl_r'])))
        SParams.append(('cntrl_Lr',float(Param_dict['cntrl_Lr'])))
    except KeyError:
        pass
    #Non-public parameters
    SParams.append(('RKbnd', 0.01))
    WParams.append(('weightrad', 5))
    WParams.append(('errbnd', 0.01))
    WParams.append(('weight0', 'Normal'))
    WParams.append(('maxrad0', 5))
    #Form dictionaries
    Std_Param = dict(SParams)
    Weight_Param = dict(WParams)
    return (Std_Param, Weight_Param)

def RasterAsc(filelocation,headlen=6):
    """Load Raster data from an ascii file, w/ or w/o header
    
    Arguments:
        - filelocation -- string
        - headlen -- header length as an int, default 6    
    Returns:
        - rastdata -- raster data as nparray 
        - header -- {header titles:header data} dictionary"""
    #load data from Maxent ascii file
    fobj = open(filelocation,'r')
    #Look for a header in the format ncols, nrows, xllcorner,
    #yllcorner, cellsize, NODATA_value
    headdata = []
    headtitles = []
    if headlen == 6:  #expected header FUTURE: make more robust
        for ii in range(0,headlen):
            line = fobj.readline()
            head = line.split()
            headtitles.append(head[0])
            if ii <= 1:
                headdata.append(int(head[1]))
            else:
                headdata.append(float(head[1]))
    fobj.close()
    #Get the actual data
    rastdata = np.loadtxt(filelocation,skiprows=headlen)
    header = dict(zip(headtitles,headdata))
    return (rastdata,header)
    
def LongLatCsv(filelocation,rowstart,rowstop,LongCol,LatCol):
    """Read a comma delimited .csv file of longitude/latitude data
    
    Arguments:
        - filelocation -- string
        - rowstart -- first row to read (first row = 1), int
        - rowstop -- last row to read, int
        - LongCol -- column with longitude data (first col = 1), int
        - LatCol -- column with latitude data, int
    Returns:
        - LongLat -- n x 2 list of longitude/latitude (x,y) tuples"""
    import csv
    
    #Read in spreadsheet, expecting long/lat data in columns
    LongLat = list(0 for ii in range(rowstop-rowstart+1)) #instantiate
    csvobj = csv.reader(open(filelocation,'rb'))
    for num, row in enumerate(csvobj):
        if (num >= rowstart-1) and (num < rowstop):
            LongLat[num-rowstart+1] = (float(row[LongCol-1]), float(row[LatCol-1]))
    return LongLat
    
def LongLatXls(filelocation,rowstart,rowstop,LongCol,LatCol,sheet=1):
    """Read an xls (2003 or earlier) file of longitude/latitude data
    
    Arguments:
        - filelocation -- string
        - rowstart -- first row to read (first row = 1), int
        - rowstop -- last row to read, int
        - LongCol -- column with longitude data (first col = 1), int
        - LatCol -- column with latitude data, int
        - sheet -- sheet index (first sheet = 1) or sheet name
    Returns:
        - LongLat -- n x 2 list of longitude/latitude (x,y) tuples"""
    if not xlrd_FLAG:
        raise ImportError('xlrd module not found: Cannot import .xls spreadsheet.')
    try:
        wb = open_workbook(filelocation)
        if type(sheet) is int:
            sheet -= 1
            pg = wb.sheet_by_index(sheet)
        elif type(sheet) is str:
            pg = wb.sheet_by_name(sheet)
        else:
            raise RuntimeError('Sheet designation must be either an int or a string')
        rowstart -= 1
        LongCol -= 1
        LatCol -= 1
        LongLat = list(0 for ii in range(rowstop-rowstart))
        for row in range(rowstart,rowstop):
            LongLat[row-rowstart] = (pg.cell(row,LongCol).value,\
                pg.cell(row,LatCol).value)
    except IOError:
        raise
    except IndexError:
        raise
    except RuntimeError:
        raise
    return LongLat

def WriteData(data, outname):
    """Write 3D data to several types of files
    
    Arguments:
        - data -- data to be written, ndarray
        - outname -- place to write the data, string"""
    for tt in range(data[0,0,:].size):
        np.savetxt(outname+'_'+str(tt)+".txt", data[:,:,tt], delimiter=',')