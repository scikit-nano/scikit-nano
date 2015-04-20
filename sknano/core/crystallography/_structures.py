# -*- coding: utf-8 -*-
"""
======================================================================
Crystal structures (:mod:`sknano.core.crystallography._structures`)
======================================================================

.. currentmodule:: sknano.core.crystallography._structures

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
#from builtins import super
#from builtins import dict
from future import standard_library
standard_library.install_aliases()
#from future.utils import with_metaclass
__docformat__ = 'restructuredtext en'

#from abc import ABCMeta, abstractproperty

#import numpy as np

#from sknano.core.math import Point, Vector
from ._base import CrystalStructure

__all__ = ['DiamondStructure', 'HexagonalClosePackedStructure',
           'CubicClosePackedStructure']


class DiamondStructure(CrystalStructure):
    """Abstract representation of diamond structure."""
    pass


class HexagonalClosePackedStructure(CrystalStructure):
    """Abstract representation of hexagonal close-packed structure."""
    pass


class CubicClosePackedStructure(CrystalStructure):
    """Abstract representation of cubic close-packed structure."""
    pass
