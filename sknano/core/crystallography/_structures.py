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

import inspect
# from abc import ABCMeta, abstractproperty

# import numpy as np

# from sknano.core.math import Point, Vector

from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from ._lattices import CrystalLattice

__all__ = ['CrystalStructure', 'DiamondStructure',
           'HexagonalClosePackedStructure',
           'CubicClosePackedStructure']


class CrystalStructure:
    """Abstract base class for crystal structures."""

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False):

        if basis is None:
            self._atoms = Atoms()
        elif basis is not None:
            self._atoms = Atoms(atoms=basis)

        self.lattice = lattice
        self.basis = basis
        self.fmtstr = "lattice={lattice!r}, basis={basis!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    @property
    def atoms(self):
        return self._atoms

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

    def __getattr__(self, name):
        if name != '_atoms':
            return getattr(self._atoms, name)

    def __setattr__(self, name, value):
        if name.startswith('_'):
            super().__setattr__(name, value)
        else:
            setattr(self._atoms, name, value)

    def __delattr__(self, name):
        if name.startswith('_'):
            super().__delattr__(name)
        else:
            delattr(self._atoms, name)

    @property
    def unit_cell(self):
        pass

    @classmethod
    def from_pymatgen(cls, *args, classmethod=None, **kwargs):
        try:
            from pymatgen import Structure
        except ImportError as e:
            print(e)
        else:
            constructor = None
            if classmethod is None:
                constructor = Structure
            else:
                constructor = getattr(Structure, classmethod, None)

            pmg_sig = inspect.signature(constructor)
            bound_sig = pmg_sig.bind(*args, **kwargs)
            print(bound_sig.signature)
            structure = constructor(*bound_sig.args, **bound_sig.kwargs)

            atoms = Atoms()
            for site in structure.sites:
                atoms.append(Atom(element=site.specie.symbol,
                                  x=site.x, y=site.y, z=site.z))
            return cls(lattice=CrystalLattice(
                       cell_matrix=structure.lattice.matrix),
                       basis=atoms)

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
