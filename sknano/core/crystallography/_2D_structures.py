# -*- coding: utf-8 -*-
"""
==============================================================================
2D crystal structures (:mod:`sknano.core.crystallography._2D_structures`)
==============================================================================

.. currentmodule:: sknano.core.crystallography._2D_structures

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.atoms import BasisAtom as Atom, BasisAtoms as Atoms
from ._2D_lattices import Crystal2DLattice
from ._base import StructureBase
from ._extras import pymatgen_structure

__all__ = ['Crystal2DStructure']


class Crystal2DStructure(StructureBase):
    """Base class for 2D crystal structures.

    .. warning:: The implementation of this class is not complete.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional
    scaling_matrix : {:class:`~python:int`, :class:`~python:list`}, optional
    structure : `Crystal3DStructure`, optional

    """
    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False,
                 scaling_matrix=None, **kwargs):
        if not isinstance(lattice, Crystal2DLattice):
            lattice = Crystal2DLattice(cell_matrix=lattice)

        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         cartesian=cartesian, scaling_matrix=scaling_matrix,
                         **kwargs)

    @classmethod
    def from_pymatgen_structure(cls, structure, scaling_matrix=None):
        atoms = Atoms()
        for site in structure.sites:
            atoms.append(Atom(element=site.specie.symbol,
                              x=site.x, y=site.y, z=site.z))
        return cls(lattice=Crystal2DLattice(
                   cell_matrix=structure.lattice.matrix), basis=atoms,
                   scaling_matrix=scaling_matrix)

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
                               classmethod='from_spacegroup',
                               scaling_matrix=scaling_matrix)
        return cls.from_pymatgen_structure(structure,
                                           scaling_matrix=scaling_matrix)
