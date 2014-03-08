# -*- coding: utf-8 -*-
"""
=============================================================
Fullerene structure tools (:mod:`sknano.nanogen._fullerenes`)
=============================================================

.. currentmodule:: sknano.nanogen._fullerenes

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

#import itertools

#from fractions import gcd
#from collections import OrderedDict

#import numpy as np

try:
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Qty = ureg.Quantity
except ImportError:
    Qty = None

#from ..chemistry import Atom

#param_units = {}

__all__ = ['Fullerene']


class Fullerene(object):
    u"""Class for creating interactive Fullerene objects.

    .. note::

       This class is under development and not functional! Attempting
       to use it raises a `RuntimeError`.

    The `fullerene isomer data
    <http://www.nanotube.msu.edu/fullerene/fullerene-isomers.html>`_
    were downloaded from the
    `The Nanotube Site <http://www.nanotube.msu.edu>`_
    at `Michigan State University <http://www.msu.edu>`_.

    Parameters
    ----------
    verbose : bool, optional
        verbose output

    Examples
    --------

    """
    def __init__(self, N=None, with_units=False, verbose=False):
        raise RuntimeError('This class is not yet implemented.')
