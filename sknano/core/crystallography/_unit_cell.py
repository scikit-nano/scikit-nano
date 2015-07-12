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

from sknano.core.atoms import BasisAtoms as Atoms

__all__ = ['UnitCell']


class UnitCell:
    """Crystal unit cell class."""

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False):

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

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

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
        self.lattice.rotate(**kwargs)
        self.basis.rotate(**kwargs)

    def todict(self):
        """Return `dict` of `UnitCell` parameters."""
        return dict(lattice=self.lattice, basis=self.basis.symbols.tolist(),
                    coords=self.basis.rs.tolist())
