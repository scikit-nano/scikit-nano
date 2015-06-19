# -*- coding: utf-8 -*-
"""
=============================================================================
Crystal lattice base class (:mod:`sknano.core.crystallography._base`)
=============================================================================

.. currentmodule:: sknano.core.crystallography._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty
from functools import total_ordering

import numpy as np

from sknano.core.math import Point
from ._unit_cell import UnitCell

__all__ = ['LatticeBase', 'ReciprocalLatticeBase', 'StructureBase']


@total_ordering
class LatticeBase:

    def __init__(self, nd=None):
        self.nd = nd
        self.offset = Point()

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def __eq__(self, other):
        if isinstance(other, type(self)):
            if self is other:
                return True
            else:
                for attr in self.__dir__():
                    if abs(getattr(self, attr) - getattr(other, attr)) > 1e-6:
                        return False
                return True
        return False

    def __lt__(self, other):
        if isinstance(other, type(self)):
            try:
                return self.cell_volume < other.cell_volume
            except AttributeError:
                return self.cell_area < other.cell_area


class ReciprocalLatticeBase(LatticeBase):

    def __getattr__(self, name):
        if name != '_direct_lattice':
            return getattr(self._direct_lattice, name)


class StructureBase:

    def __init__(self, lattice=None, basis=None, coords=None,
                 cartesian=False, scaling_matrix=None, **kwargs):

        self.unit_cell = UnitCell(lattice=lattice, basis=basis,
                                  coords=coords, cartesian=cartesian)

        self.scaling_matrix = scaling_matrix
        if scaling_matrix is not None and \
                isinstance(scaling_matrix,
                           (int, float, tuple, list, np.ndarray)):
            if isinstance(scaling_matrix, (int, float)):
                scaling_matrix = 3 * [scaling_matrix]
            self.scaling_matrix = scaling_matrix

        self.fmtstr = "{lattice!r}, {basis!r}, {coords!r}, cartesian=True"
        super().__init__(**kwargs)

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def __getattr__(self, name):
        try:
            return getattr(self.unit_cell, name)
        except AttributeError:
            return super().__getattr__(name)

    # @property
    # def atoms(self):
    #     return self._atoms

    # @atoms.setter
    # def atoms(self, value):
    #     if not isinstance(value, Atoms):
    #         raise TypeError('Expected an `Atoms` object')
    #     self._atoms = value

    @property
    def fmtstr(self):
        return self._fmtstr

    @fmtstr.setter
    def fmtstr(self, value):
        self._fmtstr = value

    def rotate(self, **kwargs):
        self.unit_cell.rotate(**kwargs)

    def todict(self):
        return self.unit_cell.todict()
