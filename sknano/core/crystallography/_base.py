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

from sknano.core import BaseClass
from sknano.core.math import Point
from ._unit_cell import UnitCell

__all__ = ['LatticeBase', 'ReciprocalLatticeBase', 'StructureBase']


@total_ordering
class LatticeBase(BaseClass):
    """Base class for crystallographic lattice objects.

    Parameters
    ----------
    nd : int
    cell_matrix : array_like
    orientation_matrix : array_like, optional

    """

    def __init__(self, nd=None, cell_matrix=None, orientation_matrix=None):
        super().__init__()

        self.nd = nd
        self.offset = Point()
        if cell_matrix is not None:
            orientation_matrix = cell_matrix.T * self.fractional_matrix

        if orientation_matrix is None:
            orientation_matrix = np.matrix(np.identity(3))

        self.orientation_matrix = orientation_matrix
        self.lattice_type = None

    def __dir__(self):
        return ['nd', 'offset', 'orientation_matrix']

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self is other or \
                all([np.allclose(getattr(self, attr), getattr(other, attr))
                     for attr in dir(self)])

    def __lt__(self, other):
        if isinstance(other, type(self)):
            try:
                return self.cell_volume < other.cell_volume
            except AttributeError:
                return self.cell_area < other.cell_area


class ReciprocalLatticeBase(LatticeBase):
    """Base class for crystallographic reciprocal lattice objects.

    Parameters
    ----------
    direct_lattice : :class:`Crystal2DLattice` or :class:`Crystal3DLattice`
    nd : int
    """
    def __init__(self, direct_lattice, nd):
        self._direct_lattice = direct_lattice
        super().__init__(
            nd=nd, cell_matrix=self._direct_lattice.cell_matrix,
            orientation_matrix=self._direct_lattice.orientation_matrix)

    def __getattr__(self, name):
        if name != '_direct_lattice':
            return getattr(self._direct_lattice, name)


@total_ordering
class StructureBase(BaseClass):
    """Base class for abstract representions of crystal structures.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional
    scaling_matrix : {:class:`~python:int`, :class:`~python:list`}, optional

    """

    def __init__(self, lattice=None, basis=None, coords=None,
                 cartesian=False, scaling_matrix=None, **kwargs):

        super().__init__(**kwargs)

        self.unit_cell = UnitCell(lattice=lattice, basis=basis,
                                  coords=coords, cartesian=cartesian)

        self.scaling_matrix = scaling_matrix
        if scaling_matrix is not None and \
                isinstance(scaling_matrix,
                           (int, float, tuple, list, np.ndarray)):
            if isinstance(scaling_matrix, (int, float)):
                scaling_matrix = 3 * [scaling_matrix]
            self.scaling_matrix = scaling_matrix

        self.fmtstr = self.unit_cell.fmtstr + \
            ", scaling_matrix={scaling_matrix!r}"

    def __dir__(self):
        return ['unit_cell', 'scaling_matrix']

    def __getattr__(self, name):
        try:
            return getattr(self.unit_cell, name)
        except AttributeError:
            return super().__getattr__(name)

    def __eq__(self, other):
        if isinstance(other, StructureBase):
            return self is other or \
                (self.unit_cell == other.unit_cell and
                 self.scaling_matrix == other.scaling_matrix)

    def __lt__(self, other):
        if isinstance(other, StructureBase):
            try:
                return ((self.scaling_matrix < other.scaling_matrix and
                         self.unit_cell <= other.unit_cell) or
                        (self.scaling_matrix <= other.scaling_matrix and
                         self.unit_cell < other.unit_cell))
            except TypeError:
                return self.unit_cell < other.unit_cell

    # @property
    # def atoms(self):
    #     return self._atoms

    # @atoms.setter
    # def atoms(self, value):
    #     if not isinstance(value, Atoms):
    #         raise TypeError('Expected an `Atoms` object')
    #     self._atoms = value

    def rotate(self, **kwargs):
        """Rotation the structure unit cell."""
        self.unit_cell.rotate(**kwargs)

    def todict(self):
        attrdict = self.unit_cell.todict()
        attrdict.update(dict(scaling_matrix=self.scaling_matrix))
        return attrdict
