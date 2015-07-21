# -*- coding: utf-8 -*-
"""
===========================================================================
Crystal unit cell class (:mod:`sknano.core.crystallography._unit_cell`)
===========================================================================

.. currentmodule:: sknano.core.crystallography._unit_cell

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering

from sknano.core import BaseClass
from sknano.core.atoms import BasisAtoms as Atoms

__all__ = ['UnitCell']


@total_ordering
class UnitCell(BaseClass):
    """Base class for abstract representations of crystallographic unit cells.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional

    """

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False):

        super().__init__()

        if basis is None:
            basis = Atoms()
        elif basis is not None:
            basis = Atoms(basis)
            basis.lattice = lattice
            if lattice is not None and coords is not None:
                for atom, pos in zip(basis, coords):
                    atom.lattice = lattice
                    if not cartesian:
                        atom.rs = pos
                    else:
                        atom.rs = lattice.cartesian_to_fractional(pos)

        self.lattice = lattice
        self.basis = basis
        self.fmtstr = "{lattice!r}, {basis!r}, {coords!r}, cartesian=False"

    def __dir__(self):
        return ['lattice', 'basis']

    def __eq__(self, other):
        return self is other or \
            all([(getattr(self, attr) == getattr(other, attr)) for attr
                 in dir(self)])

    def __lt__(self, other):
        return (self.lattice < other.lattice and self.basis <= other.basis) \
            or (self.lattice <= other.lattice and self.basis < other.basis)

    def __getattr__(self, name):
        try:
            return getattr(self.lattice, name)
        except AttributeError:
            return super().__getattr__(name)

    def __iter__(self):
        return iter(self.basis)

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

    def rotate(self, **kwargs):
        """Rotate lattice vectors and basis."""
        self.lattice.rotate(**kwargs)
        self.basis.rotate(**kwargs)

    def todict(self):
        """Return `dict` of `UnitCell` parameters."""
        return dict(lattice=self.lattice, basis=self.basis.symbols.tolist(),
                    coords=self.basis.rs.tolist())
