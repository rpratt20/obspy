# -*- coding: utf-8 -*-
"""
obspy.io.dirtree.sds - read support for SeisComP Data Structure
===============================================================
This module provides read support for data stored locally in a SeisComP Data
Structure (SDS) directory structure.

The directory and file layout of SDS is defined as:

    <SDSdir>/YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY

These fields are defined by SDS as follows:

    SDSdir :  arbitrary base directory
    YEAR   :  4 digit year
    NET    :  Network code/identifier, up to 8 characters, no spaces
    STA    :  Station code/identifier, up to 8 characters, no spaces
    CHAN   :  Channel code/identifier, up to 8 characters, no spaces
    TYPE   :  1 characters indicating the data type, recommended types are:
               'D' - Waveform data
               'E' - Detection data
               'L' - Log data
               'T' - Timing data
               'C' - Calibration data
               'R' - Response data
               'O' - Opaque data
    LOC    :  Location identifier, up to 8 characters, no spaces
    DAY    :  3 digit day of year, padded with zeros

See https://www.seiscomp3.org/wiki/doc/applications/slarchive/SDS.

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import glob
import os
from datetime import timedelta

from obspy import Stream, read

SDS_FMTSTR = os.path.join(
    "{year}", "{network}", "{station}", "{channel}.{type}",
    "{network}.{station}.{location}.{channel}.{type}.{year}.{doy:03d}")


def read_from_SDS(sds_root, seed_id, starttime, endtime, sds_type="D",
                  format="MSEED", cleanup_merge=True):
    """
    Read data from a local SeisComP Data Structure (SDS) directory tree.

    :type seed_id: str
    :param seed_id: SEED ID of data that should be read
        (e.g. 'IU.ANMO.00.LHZ'). Wildcards '*' and '?' are supported
        (e.g. 'IU.*.*.HH?').
    :type starttime: :class:`~obspy.core.utcdatetime.UTCDateTime`
    :param starttime: Start of requested time window.
    :type endtime: :class:`~obspy.core.utcdatetime.UTCDateTime`
    :param endtime: End of requested time window.
    :type sds_type: str
    :param sds_type: SDS data type identifier, one single character. Types
        recommended by the SDS definition are: 'D' - Waveform data,
        'E' - Detection data, 'L' - Log data, 'T' - Timing data,
        'C' - Calibration data, 'R' - Response data, 'O' - Opaque data. Can
        also be wildcarded by setting to `?` or `*`.
    :type format: str
    :param format: File format the data is stored in, see
        :func:`~obspy.core.stream.read()` for a list of file formats supported
        by ObsPy. Usually, SDS archives are stored in "MSEED" format. Can be
        set to `None` for file format autodetection (slowing down the reading).
    :type cleanup_merge: bool
    :param cleanup_merge: Whether to perform a cleanup merge on the resulting
        stream to merge seamless traces originating from different files
        (:meth:`st.merge(-1) <obspy.core.stream.Stream.merge>`)
    :rtype: :class:`~obspy.core.stream.Stream`
    """
    if starttime > endtime:
        msg = ("'endtime' must be after 'starttime'.")
        raise ValueError(msg)

    # SDS has data sometimes in adjacent days, so also try to read the
    # requested data from those files. Usually this is only a few seconds of
    # data after midnight, but for now we play safe here to catch all requested
    # data (and with MiniSEED - the usual SDS file format - we can use
    # starttime/endtime kwargs anyway to read only desired parts).
    year_doy = set()
    t = starttime - timedelta(days=1.01)
    t_max = endtime + timedelta(days=1.01)
    # make a list of year/doy combinations that covers the whole requested
    # time window (plus day before and day after)
    while t < t_max:
        year_doy.add((t.year, t.julday))
        t += timedelta(days=1)

    st = Stream()
    net, sta, loc, cha = seed_id.split(".")
    for year, doy in year_doy:
        filename = SDS_FMTSTR.format(
            network=net, station=sta, location=loc, channel=cha, year=year,
            doy=doy, type=sds_type)
        full_path = os.path.join(sds_root, filename)
        # check if any files match at all, read() raises an exception if not
        if not glob.glob(full_path):
            continue
        st += read(full_path, format=format, starttime=starttime,
                   endtime=endtime)

    st.trim(starttime, endtime)
    if cleanup_merge:
        st.merge(-1)
    return st


if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)
