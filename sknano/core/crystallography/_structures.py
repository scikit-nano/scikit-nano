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
# import numbers
# from abc import ABCMeta, abstractproperty

import numpy as np

# from sknano.core.math import Point, Vector

from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from ._lattices import CrystalLattice

__all__ = ['CrystalStructure', 'DiamondStructure',
           'HexagonalClosePackedStructure',
           'CubicClosePackedStructure', 'FCCStructure', 'Gold',
           'Copper', 'pymatgen_structure']


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
    def from_pymatgen_structure(cls, structure):
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


def pymatgen_structure(*args, classmethod=None, **kwargs):
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
        return constructor(*bound_sig.args, **bound_sig.kwargs)

        # atoms = Atoms()
        # for site in structure.sites:
        #     atoms.append(Atom(element=site.specie.symbol,
        #                       x=site.x, y=site.y, z=site.z))
        # return CrystalStructure(lattice=CrystalLattice(
        #                         cell_matrix=structure.lattice.matrix),
        #                         basis=atoms)


class DiamondStructure(CrystalStructure):
    """Abstract representation of diamond structure."""
    def __init__(self, a=3.567, scaling_matrix=None):
        diamond = \
            pymatgen_structure(227, CrystalLattice.cubic(a).cell_matrix, ["C"],
                               [[0, 0, 0]], classmethod='from_spacegroup')
        if scaling_matrix is not None and \
                isinstance(scaling_matrix,
                           (int, float, tuple, list, np.ndarray)):
            diamond.make_supercell(scaling_matrix)
        diamond = CrystalStructure.from_pymatgen_structure(diamond)
        super().__init__(lattice=diamond.lattice, basis=diamond.basis)


class HexagonalClosePackedStructure(CrystalStructure):
    """Abstract representation of hexagonal close-packed structure."""
    pass


class CubicClosePackedStructure(CrystalStructure):
    """Abstract representation of cubic close-packed structure."""
    pass


class FCCStructure(CrystalStructure):
    def __init__(self, lattice, basis):
        super().__init__(lattice=lattice, basis=basis)

    @classmethod
    def from_spacegroup(cls, sg, a, basis, coords, scaling_matrix=None):
        if not isinstance(basis, list):
            basis = [basis]
        if len(basis) != len(coords):
            if isinstance(coords, list) and len(coords) != 0 and \
                    isinstance(coords[0], (int, float)):
                coords = [coords]
                cls.from_spacegroup(sg, a, basis, coords,
                                    scaling_matrix=scaling_matrix)

        lattice = CrystalLattice.cubic(a)
        fcc = \
            pymatgen_structure(sg, lattice.cell_matrix, basis, coords,
                               classmethod='from_spacegroup')
        if scaling_matrix is not None and \
                isinstance(scaling_matrix,
                           (int, float, tuple, list, np.ndarray)):
            fcc.make_supercell(scaling_matrix)
        fcc = CrystalStructure.from_pymatgen_structure(fcc)
        return cls(lattice=fcc.lattice, basis=fcc.basis)


class Gold(FCCStructure):
    def __init__(self, a=4.078, scaling_matrix=None):

        gold = FCCStructure.from_spacegroup(225, a, ["Au"], [[0, 0, 0]],
                                            scaling_matrix=scaling_matrix)
        super().__init__(gold.lattice, gold.basis)


class Copper(FCCStructure):
    def __init__(self, a=3.615, scaling_matrix=None):
        copper = FCCStructure.from_spacegroup(225, a, ["Cu"], [[0, 0, 0]],
                                              scaling_matrix=scaling_matrix)
        super().__init__(copper.lattice, copper.basis)
