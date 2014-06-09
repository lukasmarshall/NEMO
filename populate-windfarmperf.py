#!/usr/bin/env python
# Load wind generation data from http://windfarmperformance.info for a
# year into the database.
#
# -*- Python -*-
# Copyright (C) 2011 Ben Elliston
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

from pylab import *
import optparse
import sys
import string
import time
import datetime
import tables

parser = optparse.OptionParser ('populate-nonsched.py')
parser.add_option("--db", type='string', default='nem.h5', help='filename')
parser.add_option("--year", type='string', help='year of data set')
parser.add_option("--compressor", type='string', default='blosc', help='PyTable compressor')
parser.add_option("--complevel", type='int', default=6, help='PyTable compression level')

opts,args = parser.parse_args ()

if not opts.year:
    parser.print_help ()
    print
    sys.exit (1)

h5file = tables.openFile(opts.db, mode = 'r+')
print h5file
try:
  h5file.createGroup(h5file.root, 'aux')
except tables.exceptions.NodeError:
  pass

try:
  h5file.createGroup(h5file.root.aux, 'windfarmperf%s' % opts.year)
except tables.exceptions.NodeError:
  print 'group windfarmperf%s already exists' % opts.year
  pass

class DispatchInterval(tables.IsDescription):
    time = tables.Time32Col(pos=0)
    duid = tables.StringCol(8,pos=1)
    power = tables.Float32Col (pos=2)

f = tables.Filters (complevel=opts.complevel, complib=opts.compressor)
table = h5file.createTable('/aux/windfarmperf%s' % opts.year, 'data', DispatchInterval, \
                               "windfarmperformance.info %s wind generation" % opts.year, filters=f)
dispatch = table.row

for month in range (12):
    f = open ('aemo_wind_%s%02d.csv' % (opts.year, month + 1), 'r')
    for count, line in enumerate (f):
	line = line.strip ()
        # "timestamp","woolnth1","captl_wf","cathrock", ...
        if count == 0:
            line = string.upper (line)
            line = string.replace (line, '"', '')
            duids = line.split (',')[1:]
        else:
            fields = line.split (',')
            timestamp = fields[0].strip ('"')
            if timestamp[0:4] != opts.year:
                print 'skipping line out of date range: ', timestamp
                continue
            t = time.mktime (time.strptime (timestamp, "%Y-%m-%d %H:%M:%S"))
            # All NEM times are UTC+10.
            t -= 60 * 60 * 10
            for duid, power in zip (duids, fields[1:]):
                if power == '':
                    # no record
                    continue
                dispatch['time'] = t
                dispatch['duid'] = duid
                dispatch['power'] = float (power)
                dispatch.append ()
    f.close ()
h5file.flush ()
