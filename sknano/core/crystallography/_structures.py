# -*- coding: utf-8 -*-
"""
===========================================================================
Crystal structure classes (:mod:`sknano.core.crystallography._structures`)
===========================================================================

.. currentmodule:: sknano.core.crystallography._structures

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
# from builtins import super
# from builtins import dict
from future import standard_library
standard_library.install_aliases()
# from future.utils import with_metaclass
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty

# import numpy as np

# from sknano.core.math import Point, Vector

from ._lattices import CrystalLattice

__all__ = ['CrystalStructure', 'DiamondStructure',
           'HexagonalClosePackedStructure',
           'CubicClosePackedStructure']


class CrystalStructure:
    """Abstract base class for crystal structures."""

    def __init__(self, lattice, basis, coords=None, cartesian=False):

        self.lattice = lattice
        self.basis = basis

        self.fmtstr = "{lattice!r}, {basis!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    @property
    def basis(self):
        """Crystal structure basis."""
        return self._basis

    @basis.setter
    def basis(self, value):
        self._basis = value

    @property
    def lattice(self):
        return self._lattice

    @lattice.setter
    def lattice(self, value):
        if not isinstance(value, CrystalLattice):
            value = CrystalLattice(cell_matrix=value)
        self._lattice = value

    @property
    def unit_cell(self):
        pass

    def todict(self):
        """Return `dict` of `CrystalStructure` parameters."""
        return dict(lattice=self.lattice, basis=self.basis)


class DiamondStructure(CrystalStructure):
    """Abstract representation of diamond structure."""
    pass


class HexagonalClosePackedStructure(CrystalStructure):
    """Abstract representation of hexagonal close-packed structure."""
    pass


class CubicClosePackedStructure(CrystalStructure):
    """Abstract representation of cubic close-packed structure."""
    pass
