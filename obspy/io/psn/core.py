# -*- coding: utf-8 -*-
"""
Public Seismic Network bindings to ObsPy core module.

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)

"""
from __future__ import (absolute_import, division, print_function, unicode_literals)

from future.builtins import *  # NOQA

from future.utils import native_str

import warnings

import struct
import binascii

import numpy as np

from obspy import Stream, Trace, UTCDateTime
from obspy.core.util import AttribDict

####################################################################

def isPSNTYPE4(filename):  
    """
    Checks whether a file is PSN format or not.

    :type filename: str
    :param filename: PSN file to be checked.
    :rtype: bool
    :return: ``True`` if a PSNVOLUME1 or a PSNTYPE4 file.
    """
    
    # file like readPSN and check for errors
    # Valid Header can be PSNVOLUME1 or PSNTYPE4
    try:
        s = struct.Struct('<10s')
        with open(filename, "rb") as fpin:
            buff = fpin.read(10)
            fpin.close
            volID = s.unpack(buff)[0]
            eventID = volID[0:8]

            # This will return false for invalid ID.
            if  (volID == 'PSNVOLUME1' or eventID == 'PSNTYPE4'):
                return True
            else:
                return False
                
    except:
        return False
    return True

########################################################################

def readPSN(filename):
    """
    Reads a PSN volume or a PSN event file and returns a stream object.

    .. warning::
        This function should NOT be called directly, it registers via the
        ObsPy :func:`~obspy.core.stream.read` function, call this instead.

    :type filename: str
    :param filename: PSN file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :returns: Stream object containing header and data.
    :current version ignores variable headers
    :current version does not parse PickInfo, seedInfo,
    :PolesandZeroes, or SensorAmptAtoD structures

    .. rubric:: Example

    >>> from obspy import read
    >>> st = read("/path/to/<filename>")
    >>> st  # doctest: +ELLIPSIS
    <obspy.core.stream.Stream object at 0x...>
    >>> print(st)  # doctest: +ELLIPSIS
    1 Trace(s) in Stream:
    54135 samples
    """
    # Initialize return data structures
    
    passthis = []

    tnumber = 0    

    # read PSN file and check for volume number of events

    with open(filename, "rb") as fpin:

        s = struct.Struct('<10sH8s')
        buff = fpin.read(20)
        volID = s.unpack(buff)[0]
        eventID = s.unpack(buff)[2]
        if volID == 'PSNVOLUME1':
            numRecords = s.unpack(buff)[1]
           
            #If volume file start reading at next byte and expect numRecords
            # of event records
            # if not volume file then re-start at beginning and read header for single event
            
            if eventID == 'PSNTYPE4':       # Volume file contains at least 1 event
                fpin.seek(12, 0)
                
        else:                               # Single event starts at beginning of file
             fpin.seek(0, 0)
             numRecords = 1
            
        while tnumber < numRecords:
            # With a Volume file the loop will start here
            # A single event file will just continue to end
            # Define event header structure, unpack and load into variables
            
            startTimestruct = struct.Struct('<H6cl')
            psnType4a = struct.Struct('<8sl')
            psnType4b = struct.Struct('<ddll3scBBddcBddd6s4s6sddhddd')
            trace = Trace()
            psn_headers = AttribDict()
            logger_info = AttribDict()
            events = AttribDict()
            
                        
            buff = fpin.read(12)    
            eventID = psnType4a.unpack(buff)[0]
            if eventID != 'PSNTYPE4':              # Verify beginning of an event record
                print( 'Bad TYPE4 header in file', test, sep=' ', end='\n')
                return False
            varHdrLength = psnType4a.unpack(buff)[1]
           
            # Read in the starttime structure

            buff = fpin.read(12)
            year, month, day, hour, minute, second, notused, \
               nanoseconds = startTimestruct.unpack(buff)
            
            # Now create obspy start time for data
            
            starttime = UTCDateTime(year, ord(month), ord(day), ord(hour), ord(minute),  \
                                   ord(second), int(nanoseconds/1000))

            # Read in the remainder of PSNTYPE4 header

            buff = fpin.read(130)        

            # Make trace.stats from header
                      
            sampling_rate = psnType4b.unpack(buff)[1]
            npts = psnType4b.unpack(buff)[2]
            station = psnType4b.unpack(buff)[15]
            channel = psnType4b.unpack(buff)[16]
            network = psnType4b.unpack(buff)[17]
            Overall_sensitivity = psnType4b.unpack(buff)[18]

            # Remainder of PSN specific header info
            
            psn_headers.timeRefType = psnType4b.unpack(buff)[4]
            psn_headers.timeRefStatus = psnType4b.unpack(buff)[5]
            psn_headers.sampleType = psnType4b.unpack(buff)[6]
            psn_headers.sampleCompression = psnType4b.unpack(buff)[7]
            psn_headers.compIncident = psnType4b.unpack(buff)[8]
            psn_headers.compAzimuth = psnType4b.unpack(buff)[9]
            psn_headers.compOrientation = psnType4b.unpack(buff)[10]
            psn_headers.sensorType = psnType4b.unpack(buff)[11]
            
            logger_info.magCorr = psnType4b.unpack(buff)[19]
            logger_info.atodBits = psnType4b.unpack(buff)[20]
            psn_headers.minimum = psnType4b.unpack(buff)[21]
            psn_headers.maximum = psnType4b.unpack(buff)[22]
            psn_headers.mean = psnType4b.unpack(buff)[23]

            trace.stats.channel = channel
            trace.stats.sampling_rate = sampling_rate
            trace.stats.starttime = starttime
            trace.stats.network = network
            trace.stats.station = station
            trace.stats.npts = npts
            trace.stats.calib = Overall_sensitivity / 100
            trace.stats.psn =  psn_headers
            trace.stats.location = u""
            trace.stats.varHdrLength = varHdrLength
            trace.stats.latitude = psnType4b.unpack(buff)[12]
            trace.stats.longitude = psnType4b.unpack(buff)[13]
            trace.stats.elevation = psnType4b.unpack(buff)[14]
            trace.stats.startTimeOffset = psnType4b.unpack(buff)[0]

            flags = psnType4b.unpack(buff)[3]        # CRC and min/max info not used

            sampleType = psnType4b.unpack(buff)[6]

            # Can be short, long int, float, or double            
            
            if sampleType == 1:
                record_dtype = np.dtype('<i4')
                
            elif sampleType == 2:
                record_dtype = np.dtype('<f4')
                
            elif sampleType == 3:
                record_dtype = np.dtype('<d8')
                
            else:
                record_dtype = np.dtype('<i2')

            # varHdr block implemented
######################################################################## 
                        
            if varHdrLength > 0:
                logger_info = AttribDict()
                event_info = AttribDict()
                                                
                trace.stats.logger = logger_info
                trace.stats.event = event_info
                                                
                vhdr = struct.Struct('<BBl')
                for l in range(0, varHdrLength):
                    buff = fpin.read(6)
                    
                    checkNumber = vhdr.unpack(buff)[0]
                    id = vhdr.unpack(buff)[1]
                    length = vhdr.unpack(buff)[2]
                    
                    if(checkNumber != 85):
                        print( 'Bad variable header in file', test, sep=' ', end='\n')
                    varHdrLength = varHdrLength - (6 + length)
                    if(id == 1):
                         trace.stats.location = fpin.read(length) 
                    elif(id == 2):
                         logger_info.sensor_info = fpin.read(length)
                    elif(id == 3):
                         event_info.comments = fpin.read(length)
                    elif(id == 6):
                         logger_info.filtering = fpin.read(length)
                    elif(id == 7):
                         logger_info.loggerID = fpin.read(length)

                    elif(id == 4):
                         if(length != 62):
                              print('Event_Info read error ', length, sep=' ',end='\n')
                         events = {}
                         details = {'latitude': 0, 'longitude': 0, 'depth': 0, 'Ms': 0, 'Mb': 0, 'Mw': 0, 'Ml': 0, 
                                 'Md': 0, 'M_oth': 0, 'other': 0,'typeCode': '','quality': '','flags': 0, 'agency': ''}
                         Timestruct = struct.Struct('<H6cl')
                         buff = fpin.read(12)
                         year, month, day, hour, minute, second, notused, \
                               nanoseconds = Timestruct.unpack(buff)
                         et = UTCDateTime(year, ord(month), ord(day), ord(hour), ord(minute),  \
                                   ord(second), int(nanoseconds/1000)) 
                         tevent = struct.Struct('<3d6h4cccH6c')
                         buff = fpin.read(length - 12)
                         details['latitude'] = tevent.unpack(buff)[0]
                         details['longitude'] = tevent.unpack(buff)[1]
                         details['depth'] = tevent.unpack(buff)[2]
                         details['Ms'] = tevent.unpack(buff)[3]
                         details['Mb'] = tevent.unpack(buff)[4]
                         details['Mw'] = tevent.unpack(buff)[5]
                         details['Ml'] = tevent.unpack(buff)[6]
                         details['Md'] = tevent.unpack(buff)[7]
                         details['M_oth'] = tevent.unpack(buff)[8]
                         details['other'] = tevent.unpack(buff)[9]
                         details['typeCode'] = tevent.unpack(buff)[10]
                         details['quality'] = tevent.unpack(buff)[11]
                         details['flags'] = tevent.unpack(buff)[12]
                         details['agency'] = tevent.unpack(buff)[13]
                         trace.stats.event[str(et)] = details
                                 
                    elif(id == 5):
                         phasedat = AttribDict()
                         trace.stats.phases = phasedat
                         pick = {}
                         if(length != 42):
                              print('Phase_pick read error ', length, sep=' ',end='\n')
                         Timestruct = struct.Struct('<H6cl')
                         buff = fpin.read(12)
                         year, month, day, hour, minute, second, notused, \
                         nanoseconds = Timestruct.unpack(buff)
                         pickt = UTCDateTime(year, ord(month), ord(day), \
                                   ord(hour), ord(minute),  \
                                   ord(second), int(nanoseconds/1000))
                         pickdat = struct.Struct('<8cHh16ch')
                         buff = fpin.read(length - 12)
                         pick['phase_name'] = pickdat.unpack(buff)[0]
                         pick['flags'] = pickdat.unpack(buff)[1]
                         pick['yloc'] = pickdat.unpack(buff)[2]
                         pick['fileName'] = pickdat.unpack(buff)[3]
                         pick['table_Depth'] = pickdat.unpack(buff)[4]
                         phasedat[str(pickt)] = pick

                    elif(id == 8 | id == 9 | id == 10):
                         international = fpin.read(length)
                         del(international)
                    
                    elif(id == 11):
                         if(length != 24):
                             print('SensorAmpAD read error ', length, sep=' ',end='\n')
                         advolts = struct.Struct('<3d')
                         buff = fpin.read(length)
                         logger_info.Sensvolts = advolts.unpack(buff)[0]
                         logger_info.gain = advolts.unpack(buff)[1]
                         logger_info.ADvolts = advolts.unpack(buff)[2]

                    elif(id == 12):             
                         paz = AttribDict()
                         poles = []
                         zeros = []
                         numpz = struct.Struct('<2H')
                         dat = struct.Struct('<dd')
                         buff = fpin.read(4)
                         paz.num_zeros = numpz.unpack(buff)[0]
                         paz.num_poles = numpz.unpack(buff)[1]
                         for zero in xrange(0, paz.num_zeros):
                              buff = fpin.read(16)
                              real = dat.unpack(buff)[0]
                              imag = dat.unpack(buff)[1]
                              zeros.append(complex(real, imag))
                         for pole in xrange(0, paz.num_poles):
                              buff = fpin.read(16)
                              real = dat.unpack(buff)[0]
                              imag = dat.unpack(buff)[1]
                              poles.append(complex(real, imag))
                         paz.zeros = zeros
                         paz.poles = poles
                         paz.constant = 100 / trace.stats.calib 
                         trace.stats.response = paz

                    elif(id == 13):
                         if(length != 8):
                             print('SeedInfo read error ', length, sep=' ',end='\n')
                         seeds = struct.Struct('<3sc3sc')
                         buff = fpin.read(length)
                         logger_info.seednet = seeds.unpack(buff)[2]
                         logger_info.seedloc = seeds.unpack(buff)[0]

                    elif(id > 13):
                         newfield = fpin.read(length)
                         print( 'Discarding unknown or new field ', length, sep=' ',end='\n')
                         del(newfield)

                    # Now at end of varHeader section of file
                    if(id == 0 & varHdrLength == 0):
                         break
 
#Pass Results of Read
                 
            indata = np.frombuffer(fpin.read(record_dtype.itemsize \
                            * npts), record_dtype)

            # Need to cast to float64 for trace object to work  Bugfix #1
            
            trace.data = indata.astype('<f8')
                                           
            # Now must discard CRC bytes and reset eventID for multiple record file 
                        
            crc = fpin.read(2)
            eventID = 'WRONG'
            passthis.append(trace)
            tnumber = tnumber + 1               # Loop count of events in file
            
            
    return Stream(traces=passthis)

