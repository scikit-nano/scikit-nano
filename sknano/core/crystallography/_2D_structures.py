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

# import inspect
# import numbers
# from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Vector

from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from ._2D_lattices import Crystal2DLattice
from ._extras import pymatgen_structure

__all__ = ['Crystal2DStructure']


class Crystal2DStructure:
    """Abstract base class for 2D crystal structures."""

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False):

        atoms = None

        if basis is None:
            atoms = Atoms()
        elif basis is not None:
            atoms = Atoms(atoms=basis)
            if coords is not None:
                for atom, pos in zip(atoms, coords):
                    if not cartesian:
                        pos = lattice.fractional_to_cartesian(pos)
                    atom.r = pos

        self._atoms = atoms
        self.lattice = lattice
        self.basis = atoms
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
        if not isinstance(value, Crystal2DLattice):
            value = Crystal2DLattice(cell_matrix=value)
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
    def from_pymatgen_structure(cls, structure):
        atoms = Atoms()
        for site in structure.sites:
            atoms.append(Atom(element=site.specie.symbol,
                              x=site.x, y=site.y, z=site.z))
        return cls(lattice=Crystal2DLattice(
                   cell_matrix=structure.lattice.matrix),
                   basis=atoms)

    @classmethod
    def from_spacegroup(cls, sg, lattice, basis, coords, scaling_matrix=None):
        if not isinstance(basis, list):
            basis = [basis]
        if len(basis) != len(coords):
            if isinstance(coords, list) and len(coords) != 0 and \
                    isinstance(coords[0], (int, float)):
                coords = [coords]
                cls.from_spacegroup(sg, lattice, basis, coords,
                                    scaling_matrix=scaling_matrix)

        structure = \
            pymatgen_structure(sg, lattice.cell_matrix, basis, coords,
                               classmethod='from_spacegroup')
        if scaling_matrix is not None and \
                isinstance(scaling_matrix,
                           (int, float, tuple, list, np.ndarray)):
            structure.make_supercell(scaling_matrix)
        structure = Crystal2DStructure.from_pymatgen_structure(structure)
        return cls(lattice=structure.lattice, basis=structure.basis)

    def todict(self):
        """Return `dict` of `Crystal2DStructure` parameters."""
        return dict(lattice=self.lattice, basis=self.basis)
