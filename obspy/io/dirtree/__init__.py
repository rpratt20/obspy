# -*- coding: utf-8 -*-
"""
obspy.io.dirtree - read support for ordered directory structures
================================================================
This module provides read support for some ordered local directory structures
(e.g. SeisComP Data Structure 'SDS'), storing data in filetypes readable by one
of ObsPy's I/O plugins (e.g. MiniSEED).

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA


if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)
