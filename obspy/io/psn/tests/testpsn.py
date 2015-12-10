#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The obspy.psn.core test suite.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import os
import unittest

from obspy import read
from obspy.core.utcdatetime import UTCDateTime
from obspy.psn.core import readpsn


class CoreTestCase(unittest.TestCase):
    """
    Test cases for psn core interface
    """
    def setUp(self):
        # directory where the test files are located
        self.path = os.path.join(os.path.dirname(__file__), 'data')

    def test_readViaObsPy(self):
        """
        Read files via obspy.core.stream.read function.
        """
        filename = os.path.join(self.path, 'testevt.psn')
        # 1
        st = read(filename)
        st.verify()
        st.sort(keys=['channel'])
        self.assertEqual(len(st), 1)
        self.assertEqual(st[0].stats.starttime,
                         UTCDateTime('2012-14-14T21:14:11.000000Z'))
        self.assertEqual(st[0].stats.endtime,
                         UTCDateTime('2012-14-14T21:16:59.990000Z'))
        self.assertEqual(len(st[0]), 21654)
        self.assertAlmostEqual(st[0].stats.sampling_rate, 6.0)
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.offset, 0.00)
        self.assertEqual(st[0].stats.Name, 'RPSD')
        self.assertEqual(st[0].stats.Network, 'PSN')
        self.assertEqual(st[0].stats.Location, 'Mitchell, SD')

    def test_readViaModule(self):
        """
        Read files via obspy.psn.core.readpsn function.
        """
        filename = os.path.join(self.path, 'testevt.psn')
        # 1
        st = readpsn(filename)
        st.verify()
        st.sort(keys=['channel'])
        self.assertEqual(len(st), 1)
        self.assertEqual(st[0].stats.starttime,
                         UTCDateTime('2012-14-14T21:14:11.000000Z'))
        self.assertEqual(st[0].stats.endtime,
                         UTCDateTime('2010-03-03T02:00:59.990000Z'))
                self.assertEqual(len(st[0]), 21654)
        self.assertAlmostEqual(st[0].stats.sampling_rate, 6.0)
        self.assertEqual(st[0].stats.channel, 'MHZ')

        self.assertEqual(st[0].stats.offset, 0.00)
        self.assertEqual(st[0].stats.Name, 'RPSD')
        self.assertEqual(st[0].stats.Network, 'PSN')
        self.assertEqual(st[0].stats.Location, 'Mitchell, SD')
        self.assertEqual(st[0].stats.latitude, '43.7109')
        self.assertEqual(st[0].stats.longitude, '-98.0083')
        self.assertEqual(st[0].stats.elevation, '400')
        self.assertEqual(st[0].stats.incident, '90')
        self.assertEqual(st[0].stats.azimuth, '90')
        self.assertEqual(st[0].stats.adbits, '16')
        self.assertEqual(st[0].stats.sensitivity, '6.35e-006')
        self.assertEqual(st[0].stats.magcorr, 'MHZ')
        self.assertEqual(st[0].stats.sensorType, 'vertical homebuilt')
        self.assertEqual(st[0].stats.depth, '4.700000')
        self.assertEqual(st[0].stats.timeref, 'WVB')
        self.assertEqual(st[0].stats.eventlatitude, '36.319000')
        self.assertEqual(st[0].stats.eventlongitude, '-96.755000')

        """
        # Open for further tests
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.channel, 'MHZ')
        self.assertEqual(st[0].stats.channel, 'MHZ')
        """
def suite():
    return unittest.makeSuite(CoreTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
