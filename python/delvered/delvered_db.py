#!/usr/bin/env python                                                                                                                                                   
import os
import sys
import numpy as np
import glob
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack, Column
from astropy.time import Time
import healpy as hp
from dlnpyutils import utils as dln, coords, db
import subprocess
import time
#from argparse import ArgumentParser
import socket
#from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
import sqlite3

def writecat2db(cat,dbfile,table='meas'):
    """ Write a catalog to the database """
    ncat = dln.size(cat)
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"S":"TEXT", "i":"INTEGER", "f":"REAL"}

    # Get the column names
    cnames = cat.dtype.names
    cdict = dict(cat.dtype.fields)
    # Create the table
    #   the primary key ROWID is automatically generated
    if len(c.execute('SELECT name from sqlite_master where type= "table" and name="'+table+'"').fetchall()) < 1:
        columns = cnames[0].lower()+' '+d2d[cdict[cnames[0]][0].kind]
        for n in cnames[1:]: columns+=', '+n.lower()+' '+d2d[cdict[n][0].kind]
        c.execute('CREATE TABLE '+table+'('+columns+')')
    # Insert statement
    columns = []
    for n in cnames: columns.append(n.lower())
    qmarks = np.repeat('?',dln.size(cnames))
    c.executemany('INSERT INTO '+table+'('+','.join(columns)+') VALUES('+','.join(qmarks)+')', list(cat))
    db.commit()
    db.close()

def createindexdb(dbfile,col='measid',table='meas',unique=True,verbose=False):
    """ Index a column in the database """
    t0 = time.time()
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    index_name = 'idx_'+col+'_'+table
    # Check if the index exists first
    c.execute('select name from sqlite_master')
    d = c.fetchall()
    for nn in d:
        if nn[0]==index_name:
            print(index_name+' already exists')
            return
    # Create the index
    if verbose: print('Indexing '+col)
    if unique:
        c.execute('CREATE UNIQUE INDEX '+index_name+' ON '+table+'('+col+')')
    else:
        c.execute('CREATE INDEX '+index_name+' ON '+table+'('+col+')')
    data = c.fetchall()
    db.close()
    if verbose: print('indexing done after '+str(time.time()-t0)+' sec')

def getdatadb(dbfile,table='meas',cols='*',rar=None,decr=None,verbose=False):
    """ Get rows from the database """
    t0 = time.time()
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"TEXT":(np.str,100), "INTEGER":np.int, "REAL":np.float}

    # Start the SELECT statement
    cmd = 'SELECT '+cols+' FROM '+table
    # RA constraints
    if rar is not None:
        if cmd.find('WHERE') == -1:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        cmd += 'ra>='+str(rar[0])+' AND ra<'+str(rar[1])
    # DEC constraints
    if decr is not None:
        if cmd.find('WHERE') == -1:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        cmd += 'dec>='+str(decr[0])+' AND dec<'+str(decr[1])

    # Execute the select command
    #print('CMD = '+cmd)
    cur.execute(cmd)
    data = cur.fetchall()

    # No results
    if len(data)==0:
        return np.array([])

    # Get table column names and data types
    cur.execute("select sql from sqlite_master where tbl_name = '"+table+"'")
    dum = cur.fetchall()
    db.close()
    head = dum[0][0]
    # 'CREATE TABLE exposure(expnum TEXT, nchips INTEGER, filter TEXT, exptime REAL, utdate TEXT, uttime TEXT, airmass REAL, wcstype TEXT)'
    lo = head.find('(')
    hi = head.find(')')
    head = head[lo+1:hi]
    cols = head.split(',')
    cols = dln.strip(cols)
    dt = []
    for c in cols:
        pair = c.split(' ')
        dt.append( (pair[0], d2d[pair[1]]) )
    dtype = np.dtype(dt)

    # Convert to nump structured array
    cat = np.zeros(len(data),dtype=dtype)
    cat[...] = data
    del(data)

    if verbose: print('got data in '+str(time.time()-t0)+' sec.')

    return cat

def createsumtable(dbfile=None,delvedir='/net/dl2/dnidever/delve/'):
    """ Create the DELVE-MC exposures and chips database"""

    expdir = delvedir+'exposures/'
    if dbfile is None:
        dbfile = delvedir+'bricks/db/delvered_summary.db'

    # Loop over list of nightly summary files

    # Get all of the nights
    dirs = glob.glob(expdir+'20??????')
    dirs = np.array(dirs)
    ndirs = dln.size(dirs)
    nights = np.zeros(dirs.size,(np.str,10))
    for i,d in enumerate(dirs): nights[i]=os.path.basename(d)
    si = np.argsort(nights)
    nights = nights[si]
    dirs = dirs[si]
    # Nightly summary files
    allnightsumfiles = dln.strjoin(dirs,'/',nights)
    allnightsumfiles = dln.strjoin(allnightsumfiles,'_summary.fits')
    gdnightsumfiles,nnightsumfiles = dln.where(dln.exists(allnightsumfiles))
    print(str(nnightsumfiles)+' nightly summary files found')
    nightsumfiles = allnightsumfiles[gdnightsumfiles]
    # Load in the data
    expcount = 0
    chcount = 0
    lexpstr = []  # lists
    lchstr = []
    for i in range(nnightsumfiles):
        print('Loading '+nightsumfiles[i])
        #st = os.stat(nightsumfiles[i])
        try:
            expstr1 = fits.getdata(nightsumfiles[i],1,memmap=False)
            nexpstr1 = dln.size(expstr1)
            lexpstr.append(expstr1)
            chstr1 = fits.getdata(nightsumfiles[i],2,memmap=False)
            nchstr1 = dln.size(chstr1)
            lchstr.append(chstr1)
        except:
            print('Problem loading '+nightsumfiles[i])
        
       # # Load the database
       # writecat2db(expstr1,dbfile,'exposure')
       # writecat2db(chstr1,dbfile,'chip')

    # Concatenate them
    expstr = dln.concatenate(lexpstr)
    chstr = dln.concatenate(lchstr)

    # Load the database
    writecat2db(expstr,dbfile,'exposure')
    writecat2db(chstr,dbfile,'chip')


    # Create the indices
    print('Indexing')
    createindexdb(dbfile,'ra','exposure',unique=False)
    createindexdb(dbfile,'dec','exposure',unique=False)
    db.analyzetable(dbfile,'exposure')
    createindexdb(dbfile,'ra','chip',unique=False)
    createindexdb(dbfile,'dec','chip',unique=False)
    db.analyzetable(dbfile,'chip')
