
#
# Copyright (C) 2009 Tobias Megies, Moritz Beyreuther, Yannik Behr
# --------------------------------------------------------------------
"""
Various Seismogram Filtering Functions
:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
"""

import warnings

import numpy as np
from scipy.fftpack import hilbert
from scipy.signal import (cheb2ord, cheby2, convolve, get_window, iirfilter,
                          remez, decimate)



# simple downsampling using scipy.signal.decimate
def downsample (olddata,oldrate,newrate):  
  ratio=oldrate/float(newrate) # Calculate ratio of sampling rates
  newdata = decimate(olddata, int(ratio), ftype='fir',zero_phase=True)
  return newdata    