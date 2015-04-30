# -*- coding: utf-8 -*-
"""
======================================================================
Lattice systems (:mod:`sknano.core.crystallography._lattices`)
======================================================================

.. currentmodule:: sknano.core.crystallography._lattices

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
from builtins import super
from builtins import dict
from future import standard_library
standard_library.install_aliases()
# from future.utils import with_metaclass
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Point, Vector

__all__ = ['BravaisLattice', 'SimpleCubicLattice',
           'BodyCenteredCubicLattice', 'FaceCenteredCubicLattice']


class BravaisLattice:
    """Class for bravais lattices."""
    pass


class SimpleCubicLattice(BravaisLattice):
    """Abstract representation of simple cubic lattice."""
    pass


class BodyCenteredCubicLattice(BravaisLattice):
    """Abstract representation of body-centered cubic lattice."""
    pass


class FaceCenteredCubicLattice(BravaisLattice):
    """Abstract representation of face-centered cubic lattice."""
    pass
