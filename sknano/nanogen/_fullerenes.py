# -*- coding: utf-8 -*-
"""
=============================================================
Fullerene structure tools (:mod:`sknano.nanogen._fullerene`)
=============================================================

.. currentmodule:: sknano.nanogen._fullerene

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

#import itertools
import warnings
warnings.filterwarnings('ignore')  # to suppress the Pint UnicodeWarning

#from fractions import gcd
#from collections import OrderedDict

#import numpy as np

try:
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Qty = ureg.Quantity
except ImportError:
    Qty = None

from ..chemistry import Atom
from ..tools.refdata import CCbond

#param_units = {}

__all__ = ['Fullerene']


class Fullerene(object):
    u"""Class for creating interactive Fullerene objects.

    Parameters
    ----------
    verbose : bool, optional
        verbose output

    Examples
    --------

    """
    def __init__(self, N=None, with_units=False, verbose=False):
        raise RuntimeError('This class is not yet implemented.')
