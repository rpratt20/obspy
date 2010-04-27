# -*- coding: utf-8 -*-
"""
obspy.xseed: Tool to convert between Dataless SEED and XML-SEED files
=====================================================================
XML-SEED was introduced by Tsuboi, Tromp and Komatitsch (2004), it is a XML
representation of Dataless SEED. This module contains converters from Dataless
SEED to XML-SEED and vice versa as well as a converter from Dataless SEED to
RESP files. The obspy.xseed module is tested against the complete ORFEUS
Dataless SEED archive, the IRIS (US) Dataless SEED archive and against ArcLink
response requests.  

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)

Allocate a Parser object and read/write
---------------------------------------
>>> from obspy.xseed import Parser

>>> sp = Parser("dataless.seed.BW_RJOB")
>>> sp.writeXSEED("dataless.seed.BW_RJOB.xml")

The lines above will convert Dataless SEED, e.g.::

    000001V 010009402.3121970,001,00:00:00.0000~2038,001,00:00:00.0000~
    2009,037,04:32:41.0000~BayernNetz~~0110032002RJOB 000003RJOB 000008
    ...

to the XML-SEED representation, e.g.::

    <?xml version='1.0' encoding='utf-8'?>
    <xseed version="1.0">
      <volume_index_control_header>
        <volume_identifier blockette="010">
          <version_of_format>2.4</version_of_format>
          <logical_record_length>12</logical_record_length>
          <beginning_time>1970-01-01T00:00:00</beginning_time>
          <end_time>2038-01-01T00:00:00</end_time>
          <volume_time>2009-02-06T04:32:41</volume_time>
          <originating_organization>BayernNetz</originating_organization>
    ...


A response file can be written in a similar manner, just replace writeXSEED
by writeRESP:

>>> sp.writeRESP(folder="BW_RJOB", zipped=False)



xseed.Parser Object
-------------------
Every SEED and XSEED object will be parsed in a xseed.Parser structure. ObsPy
currently only supports dataless SEED/XSEED values but this will hopefully
change soon.

SEED volumes have four different volume types:

* Volume Index Control Headers
* Abbreviation Dictionary Control Headers
* Station Control Headers
* Time Span Control Headers (currently not supported by ObsPy. Some dummy
  headers will be written in case they are needed by SEED/XSEED conventions.) 

After parsing a SEED/XSEED file the Blockette objects for each volume will be
stored in Parser.volume, Parser.abbreviations and Parser.stations. Each item is
a list that contains all Blockettes and Parser.stations is a list that contain
a list for each station. 
"""

# needs to stay above import statements
DEFAULT_XSEED_VERSION = '1.1'

from obspy.core.util import _getVersionString
from obspy.xseed.parser import Parser


__version__ = _getVersionString("obspy.xseed")
