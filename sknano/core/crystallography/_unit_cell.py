# -*- coding: utf-8 -*-
"""
===========================================================================
Crystal unit cell classes (:mod:`sknano.core.crystallography._unit_cells`)
===========================================================================

.. currentmodule:: sknano.core.crystallography._unit_cells

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
# from builtins import super
# from builtins import dict
from future import standard_library
standard_library.install_aliases()
# from future.utils import with_metaclass

__docformat__ = 'restructuredtext en'

import copy
# import numbers
# from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Vector
from sknano.core.atoms import BasisAtoms as Atoms

__all__ = ['UnitCell']


class UnitCell:
    """Abstract base class for crystal unit cell."""

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False):

        if basis is None:
            basis = Atoms()
        elif basis is not None:
            basis = Atoms(basis, lattice=lattice)
            if coords is not None:
                for atom, pos in zip(basis, coords):
                    atom.lattice = lattice
                    if not cartesian:
                        atom.rs = pos
                    else:
                        atom.rs = lattice.cartesian_to_fractional(pos)

        self.lattice = lattice
        self.basis = basis
        self.fmtstr = "{lattice!r}, {basis!r}, {coords!r}, cartesian=False"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    # @property
    # def basis(self):
    #     return self._basis

    # @basis.setter
    # def basis(self, value):
    #     lattice = self.lattice
    #     coords = self.coords
    #     if value is None:
    #         value = Atoms()
    #     elif value is not None:
    #         value = Atoms(value, lattice=lattice)
    #         if coords is not None:
    #             for atom, pos in zip(basis, coords):
    #                 atom.lattice = lattice
    #                 if not cartesian:
    #                     atom.rs = pos
    #                 else:
    #                     atom.rs = lattice.cartesian_to_fractional(pos)

    # @property
    # def lattice(self):
    #     return self._lattice

    # @lattice.setter
    # def lattice(self, value):
    #     self._lattice = value

    def __getattr__(self, name):
        try:
            return getattr(self.lattice, name)
        except AttributeError:
            return super().__getattr__(name)

        # if name != '_lattice':
        #     return getattr(self._lattice, name)

    # def __setattr__(self, name, value):
    #     if name.startswith('_'):
    #         super().__setattr__(name, value)
    #     else:
    #         setattr(self._lattice, name, value)

    def __iter__(self):
        return iter(self.basis)

    def rotate(self, **kwargs):
        self.lattice.rotate(**kwargs)
        self.basis.rotate(**kwargs)

    def todict(self):
        """Return `dict` of `CrystalStructure` parameters."""
        return dict(lattice=self.lattice, basis=self.basis.symbols.tolist(),
                    coords=self.basis.rs.tolist())
