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

# import numbers
# from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Vector

# from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.core.atoms import BasisAtom as Atom, BasisAtoms as Atoms
from ._3D_lattices import CrystalLattice
from ._base import StructureBase
from ._extras import pymatgen_structure

__all__ = ['CrystalStructure', 'Crystal3DStructure', 'DiamondStructure',
           'HexagonalClosePackedStructure', 'CubicClosePackedStructure',
           'BCCStructure', 'FCCStructure', 'HexagonalStructure',
           'Gold', 'Copper', 'AlphaQuartz']


class Crystal3DStructure(StructureBase):
    """Abstract base class for crystal structures."""

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False):
        if not isinstance(lattice, CrystalLattice):
            lattice = CrystalLattice(cell_matrix=lattice)

        super().__init__(lattice, basis, coords, cartesian)

    @classmethod
    def from_pymatgen_structure(cls, structure):
        atoms = Atoms()
        for site in structure.sites:
            atoms.append(Atom(element=site.specie.symbol,
                              x=site.x, y=site.y, z=site.z))
        return cls(lattice=CrystalLattice(
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
        structure = CrystalStructure.from_pymatgen_structure(structure)
        return cls(lattice=structure.lattice, basis=structure.basis)

    def todict(self):
        """Return `dict` of `CrystalStructure` parameters."""
        return dict(lattice=self.lattice, basis=self.basis)

CrystalStructure = Crystal3DStructure


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


class HexagonalStructure(CrystalStructure):
    @classmethod
    def from_spacegroup(cls, sg, a, c, basis, coords, scaling_matrix=None):
        lattice = CrystalLattice.hexagonal(a, c)
        # lattice = CrystalLattice(a=a, b=a, c=c, alpha=90, beta=90, gamma=60)
        return CrystalStructure.from_spacegroup(sg, lattice, basis, coords,
                                                scaling_matrix=scaling_matrix)


class BCCStructure(CrystalStructure):
    pass


class FCCStructure(CrystalStructure):
    @classmethod
    def from_spacegroup(cls, sg, a, basis, coords, scaling_matrix=None):
        lattice = CrystalLattice.cubic(a)
        return CrystalStructure.from_spacegroup(sg, lattice, basis, coords,
                                                scaling_matrix=scaling_matrix)


class AlphaQuartz(HexagonalStructure):
    def __init__(self, a=4.916, c=5.405, scaling_matrix=None):
        lattice = CrystalLattice.hexagonal(a, c)
        basis = 3 * ["Si"]
        basis.extend(6 * ["O"])
        basis = Atoms(atoms=[Atom(element=s) for s in basis])
        coords = [[0.4697, 0.0000, 0.0000],
                  [0.0000, 0.4697, 0.6667],
                  [0.5305, 0.5303, 0.3333],
                  [0.4133, 0.2672, 0.1188],
                  [0.2672, 0.4133, 0.5479],
                  [0.7328, 0.1461, 0.7855],
                  [0.5867, 0.8539, 0.2145],
                  [0.8539, 0.5867, 0.4521],
                  [0.1461, 0.7328, 0.8812]]
        structure = pymatgen_structure(lattice.cell_matrix,
                                       basis.symbols,
                                       coords)
        if scaling_matrix is not None and \
                isinstance(scaling_matrix,
                           (int, float, tuple, list, np.ndarray)):
            structure.make_supercell(scaling_matrix)
        structure = CrystalStructure.from_pymatgen_structure(structure)
        super().__init__(structure.lattice, structure.basis)
        # super().__init__(lattice, basis, coords)

        # quartz = \
        #     HexagonalStructure.from_spacegroup(154, a, c, ["Si", "O"],
        #                                        [[0.4697, 0.0000, 0.0000],
        #                                         [0.4135, 0.2669, 0.1191]],
        #                                        scaling_matrix=scaling_matrix)
        # super().__init__(quartz.lattice, quartz.basis)


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
