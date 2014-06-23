# -*- Python -*-
# Copyright (C) 2011, 2013 Ben Elliston
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# If you find a bug or implement an enhancement, please send a patch
# to the author in the form of a unified context diff (diff -u) to
# <b.elliston@student.unsw.edu.au>.

import os
import bz2
import math
import sys
import argparse
import datetime
from latlong import LatLong

# PyEphem, from http://rhodesmill.org/pyephem/
# PyEphem provides scientific-grade astronomical computations
import ephem

# From Paul Gilman <Solar.Advisor.Support@nrel.gov>:
# The first list shows the data columns SAM reads from the weather file:

# Dry bulb temperature
# Dew point temperature
# Wet bulb temperature
# Percent relative humidity
# Wind velocity
# Wind direction
# Atmospheric pressure
# Global horizontal radiation (not interpolated)
# Direct normal radiation (not interpolated)
# Latitude
# Longitude
# Site elevation
# Hour of the day

### Sample row from the BoM weather data:
### hm, 48027,2009,01,01,00,00, 22.3,N, 13.1,N,  3.4,N, 29,N,   7.8,N,  26.9,N,  1.5,N,220,N,  2.1,N,1005.6,N, 975.6,N, 1,#


def verbose(s):
    if args.verbose:
        print >>sys.stderr, s


def warn(s):
    if args.verbose:
        print >>sys.stderr, 'warning:',
        print >>sys.stderr, s


def verify(items):
    # Verify that the line is valid.
    if items[0] != 'hm':
        warn('non-hm record')

    st = items[1].strip().lstrip('0')
    if st != stnumber:
        print '%s is a foreign station number' % st


def tmy3_preamble(f):
    # eg. 722287,"ANNISTON METROPOLITAN AP",AL,-6.0,33.583,-85.850,186
    print >>f, '%s in %s,\"%s\",%s,%.1f,%.3f,%.3f,%d' % \
        (stnumber, stname, args.year, ststate[0:2], args.tz, locn.lat, locn.lon, elevation)
    print >>f, 'Date (MM/DD/YYYY),Time (HH:MM),ETR (W/m^2),ETRN (W/m^2),GHI (W/m^2),GHI source,GHI uncert (%),DNI (W/m^2),DNI source,DNI uncert (%),DHI (W/m^2),DHI source,DHI uncert (%),GH illum (lx),GH illum source,Global illum uncert (%),DN illum (lx),DN illum source,DN illum uncert (%),DH illum (lx),DH illum source,DH illum uncert (%),Zenith lum (cd/m^2),Zenith lum source,Zenith lum uncert (%),TotCld (tenths),TotCld source,TotCld uncert (code),OpqCld (tenths),OpqCld source,OpqCld uncert (code),Dry-bulb (C),Dry-bulb source,Dry-bulb uncert (code),Dew-point (C),Dew-point source,Dew-point uncert (code),RHum (%),RHum source,RHum uncert (code),Pressure (mbar),Pressure source,Pressure uncert (code),Wdir (degrees),Wdir source,Wdir uncert (code),Wspd (m/s),Wspd source,Wspd uncert (code),Hvis (m),Hvis source,Hvis uncert (code),CeilHgt (m),CeilHgt source,CeilHgt uncert (code),Pwat (cm),Pwat source,Pwat uncert (code),AOD (unitless),AOD source,AOD uncert (code),Alb (unitless),Alb source,Alb uncert (code),Lprecip depth (mm),Lprecip quantity (hr),Lprecip source,Lprecip uncert (code)'


def epw_preamble(f):
    print >>f, 'LOCATION,%s (%s) in %s,%s,AUS,BoM,%s,%.2f,%.2f,%.1f,%.1f' % \
        (stname, stnumber, args.year, ststate, stnumber, locn.lat, locn.lon, opts.args, elevation)

    print >>f, 'DESIGN CONDITIONS,0'
    print >>f, 'TYPICAL/EXTREME PERIODS,,'
    print >>f, 'GROUND TEMPERATURES,,,,,,'
    print >>f, 'HOLIDAYS/DAYLIGHT SAVINGS,No,0,0,0'
    print >>f, 'COMMENTS 1,Generated by weather-maker.py from Bureau of Meteorology solar and weather data (%d)' % opts.year
    print >>f, 'COMMENTS 2,Please report bugs in weather-maker.py to b.elliston@student.unsw.edu.au'
    print >>f, 'DATA PERIODS,1,1,Data,Sunday,1/ 1,12/31'


def tmy3_record(f, record):
    t = datetime.datetime(args.year, 1, 1)
    t += datetime.timedelta(hours=record['hour'])

    line = '%02d/%02d/%d,%02d:50,-9900,-9900,%d,1,5,%d,1,5,-9900,1,0,-9900,1,0,-9900,1,0,-9900,1,0,-9900,1,0,-9900,?,9,-9900,?,9,%.1f,A,7,%.1f,A,7,%.1f,A,7,%d,A,7,%d,A,7,%.1f,A,7,-9900,?,9,-9900,?,9,-9900,?,9,-9900,?,9,-9900,?,9,-9900,-9900,?,9' \
        % (t.month, t.day, t.year, t.hour + 1, record['ghi'], record['dni'],
           record['dry-bulb'], record['dew-point'], record['rel-humidity'],
           record['atm-pressure'] / 100, record['wind-direction'], record['wind-speed'])
    print >>f, line


def epw_record(f, record):
    t = datetime.datetime(args.year, 1, 1)
    t += datetime.timedelta(hours=record['hour'])

    line = '%d,%d,%d,%d,50,_______________________________________,%.1f,%.1f,%d,%d,9999,9999,9999,%d,%d,%d,999999,999999,999999,999999,%d,%.1f,99,99,9999,99999,9,999999999,99999,0.999,999,99,999,0,99' \
        % (t.year, t.month, t.day, t.hour + 1, record['dry-bulb'], record['dew-point'],
           record['rel-humidity'], record['atm-pressure'], record['ghi'], record['dni'], record['dhi'],
           record['wind-direction'], record['wind-speed'])
    print >>f, line


# Return the GHI and DNI for a given location and time.
def irradiances(locn, hour):
    x, y = locn.xy()
    # Compute a solar data filename from the hour
    # Use 2010 as the reference year, as it was not a leap year.
    hours = datetime.timedelta(hours=hour)
    tzoffset = datetime.timedelta(hours=args.tz)
    hr = datetime.datetime(args.year, 1, 1) + hours - tzoffset
    if hr.month == 2 and hr.day == 29:
        # skip Feb 29 on leap years
        hr += datetime.timedelta(days=1)
    filename = hr.strftime(args.grids + '/HOURLY_GHI/%d/' % args.year + hr.strftime('solar_ghi_%Y%m%d_%HUT.txt'))
    try:
        f = bz2.BZ2File(filename + '.bz2', 'r')
        line = f.readlines()[x + 6]
        f.close()
        ghr = int(line.split()[y])
    except IOError:
        try:
            f = open(filename, 'r')
            line = f.readlines()[x + 6]
            f.close()
            ghr = int(line.split()[y])
        except IOError:
            # print 'missing', filename
            ghr = 0

    filename = hr.strftime(args.grids + '/HOURLY_DNI/%d/' % args.year + hr.strftime('solar_dni_%Y%m%d_%HUT.txt'))
    try:
        f = bz2.BZ2File(filename + '.bz2', 'r')
        line = f.readlines()[x + 6]
        f.close()
        dnr = int(line.split()[y])
    except IOError:
        try:
            f = open(filename, 'r')
            line = f.readlines()[x + 6]
            f.close()
            dnr = int(line.split()[y])
        except IOError:
            # print 'missing', filename
            dnr = 0

    if ghr == -999:
        ghr = 0
    if dnr == -999:
        dnr = 0

    # Compute direct horizontal irradiance:
    # DHI = GHI - DNI cos (zenith)
    observer.date = hr + datetime.timedelta(minutes=50)
    sun.compute(observer)
    zenith = (math.pi / 2.) - sun.alt
    dhr = ghr - dnr * math.cos(zenith)
    if dhr < -10:
        # Don't worry about diffuse levels below 10 W/m2.
        warn('negative diffuse horizontal irradiance: %d' % dhr)
        dhr = 0
    return ghr, dnr, dhr


# Read station details file.
def station_details():
    global stnumber
    global stname
    global ststate

    line = [line for line in open(args.hm_details, 'r') if 'st,' + args.st in line][0]
    st = line[0:2]
    stnumber = line[3:9].strip().lstrip('0')
    stname = line[15:55].strip()
    ststate = line[107:110]
    verbose('Processing station number %s (%s)' % (stnumber, stname))

    latitude = float(line[72:80])
    longitude = float(line[81:90])
    locn = LatLong((latitude, longitude))
    altitude = int(float(line[111:117]))
    wflags = line[153:156]
    sflags = line[157:160]
    iflags = line[161:164]
    if int(wflags) or int(sflags) or int(iflags):
        warn('%% wrong = %s, %% suspect = %s, %% inconsistent = %s'
             % (wflags, sflags, iflags))

    return(locn, altitude)

parser = argparse.ArgumentParser(description='Bug reports to: b.elliston@student.unsw.edu.au')
parser.add_argument('--version', action='version', version='1.1')
parser.add_argument("--grids", type=str, help='top of gridded data tree', required=True)
parser.add_argument("-y", "--year", type=int, help='year to generate', required=True)
parser.add_argument("--st", type=str, help='BoM station code (required)', required=True)
parser.add_argument("--hm-data", type=str, help='BoM station data file', required=True)
parser.add_argument("--hm-details", type=str, help='BoM station details file', required=True)
parser.add_argument("--tz", type=float, default=10.0, help='Time zone [default +10]')
parser.add_argument("-o", "--out", type=str, help='output filename', required=True)
parser.add_argument("--format", "--format", default="epw", help="output format: EPW [default], TMY3", metavar="FORMAT")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="verbose run output")
args = parser.parse_args()

# Check that the grid directory exists
if not os.path.isdir(args.grids):
    print >>sys.stderr, 'error: %s is not a directory' % args.grids
    sys.exit(1)

infile = open(args.hm_data, 'r')
outfile = open(args.out, 'wb')

locn, elevation = station_details()
sun = ephem.Sun()
observer = ephem.Observer()
observer.elevation = elevation
observer.lat = str(locn.lat)
observer.long = str(locn.lon)

if args.format.lower() == 'tmy3':
    verbose('Generating a TMY3 file')
    tmy3_preamble(outfile)
elif args.format.lower() == 'epw':
    verbose('Generating an EPW file')
    epw_preamble(outfile)
else:
    raise ValueError("unknown format %s" % args.format)

i = 0
for line in infile:
    if len(line) == 1:
        # Skip weird ^Z lines.
        continue
    data = line.split(',')
    if data[1] == 'Station Number':
        # Skip this line; it is the header.
        continue
    if data[2] != str(args.year):
        # Skip years that are not of interest.
        continue
    if data[3] == '02' and data[4] == '29':
        warn('skipping Feb 29')
        i += 1
        continue

    # Generate pedantic warnings.
    verify(data)

    record = {}
    record['hour'] = i
    try:
        record['dry-bulb'] = float(data[7])
    except ValueError:
        record['dry-bulb'] = 99.9
    try:
        record['wet-bulb'] = float(data[8])
    except ValueError:
        record['wet-bulb'] = 99.9
    try:
        record['dew-point'] = float(data[9])
    except ValueError:
        record['dew-point'] = 99.9
    try:
        record['rel-humidity'] = float(data[10])
    except ValueError:
        record['rel-humidity'] = 999.
    try:
        record['wind-speed'] = float(data[11])
    except ValueError:
        record['wind-speed'] = 999.
    try:
        record['wind-direction'] = int(data[12])
    except ValueError:
        record['wind-direction'] = 999
    try:
        record['atm-pressure'] = int(float(data[15]) * 100)
    except ValueError:
        record['atm-pressure'] = 999999.

    record['ghi'], record['dni'], record['dhi'] = irradiances(locn, i)
    i += 1

    if args.format.lower() == 'tmy3':
        tmy3_record(outfile, record)
    elif args.format.lower() == 'epw':
        epw_record(outfile, record)

infile.close()
outfile.close()
