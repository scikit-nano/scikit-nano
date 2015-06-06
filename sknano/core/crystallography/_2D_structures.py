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

# from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.core.atoms import BasisAtom as Atom, BasisAtoms as Atoms
from ._2D_lattices import Crystal2DLattice
from ._base import StructureBase
from ._extras import pymatgen_structure

__all__ = ['Crystal2DStructure']


class Crystal2DStructure(StructureBase):
    """Abstract base class for 2D crystal structures."""

    def __init__(self, lattice=None, basis=None, coords=None,
                 cartesian=False):
        if not isinstance(lattice, Crystal2DLattice):
            lattice = Crystal2DLattice(cell_matrix=lattice)

        super().__init__(lattice, basis, coords, cartesian)

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
