# -*- coding: utf-8 -*-
"""
===========================================================================
Crystal structure classes (:mod:`sknano.core.crystallography._structures`)
===========================================================================

.. currentmodule:: sknano.core.crystallography._structures

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.atoms import BasisAtom as Atom, BasisAtoms as Atoms
from sknano.core.refdata import lattice_parameters as lattparams
from ._3D_lattices import Crystal3DLattice
from ._base import StructureBase
from ._extras import pymatgen_structure

__all__ = ['CrystalStructure', 'Crystal3DStructure',
           'CaesiumChlorideStructure', 'CsClStructure',
           'DiamondStructure',
           'RocksaltStructure', 'RockSaltStructure', 'NaClStructure',
           'SphaleriteStructure', 'ZincblendeStructure', 'ZincBlendeStructure',
           'BCCStructure',
           'FCCStructure', 'Copper', 'Gold',
           'CubicClosePackedStructure', 'CCPStructure',
           'HexagonalClosePackedStructure', 'HCPStructure',
           'HexagonalStructure', 'AlphaQuartz', 'MoS2']


class Crystal3DStructure(StructureBase):
    """Abstract base class for crystal structures."""

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False,
                 scaling_matrix=None, structure=None, **kwargs):

        if structure is not None:
            lattice = structure.lattice
            basis = structure.basis
            scaling_matrix = structure.scaling_matrix

        if not isinstance(lattice, Crystal3DLattice):
            lattice = Crystal3DLattice(cell_matrix=lattice)

        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         cartesian=cartesian, scaling_matrix=scaling_matrix,
                         **kwargs)

    @classmethod
    def from_pymatgen_structure(cls, structure, scaling_matrix=None):
        atoms = Atoms()
        for site in structure.sites:
            atoms.append(Atom(site.specie.symbol,
                              x=site.x, y=site.y, z=site.z))
        return cls(lattice=Crystal3DLattice(
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

CrystalStructure = Crystal3DStructure


class MoS2(Crystal3DStructure):
    """Molybdenum disulphide structure class."""
    def __init__(self, a=lattparams['molybdenum_disulphide']['a'],
                 c=lattparams['molybdenum_disulphide']['c'],
                 basis=['Mo', 'S'], scaling_matrix=None, **kwargs):

        molybdenum_disulphide = \
            pymatgen_structure(194,
                               Crystal3DLattice.hexagonal(a, c).cell_matrix,
                               basis, [[1/3, 2/3, 1/4], [1/3, 2/3, 0.621]],
                               classmethod='from_spacegroup',
                               scaling_matrix=scaling_matrix)
        molybdenum_disulphide = \
            Crystal3DStructure.from_pymatgen_structure(
                molybdenum_disulphide, scaling_matrix=scaling_matrix)

        super().__init__(structure=molybdenum_disulphide, **kwargs)


class CaesiumChlorideStructure(Crystal3DStructure):
    """Abstract representation of caesium chloride structure."""
    def __init__(self, a=lattparams['caesium_chloride'], basis=['Cs', 'Cl'],
                 scaling_matrix=None, **kwargs):
        caesium_chloride = \
            pymatgen_structure(221, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0], [0.5, 0.5, 0.5]],
                               classmethod='from_spacegroup',
                               scaling_matrix=scaling_matrix)
        caesium_chloride = \
            Crystal3DStructure.from_pymatgen_structure(
                caesium_chloride, scaling_matrix=scaling_matrix)

        super().__init__(structure=caesium_chloride, **kwargs)

CsClStructure = CaesiumChlorideStructure


class DiamondStructure(Crystal3DStructure):
    """Abstract representation of diamond structure."""
    def __init__(self, a=lattparams['diamond'], basis=['C'],
                 scaling_matrix=None, **kwargs):
        diamond = \
            pymatgen_structure(227, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0]],
                               classmethod='from_spacegroup',
                               scaling_matrix=scaling_matrix)
        diamond = \
            Crystal3DStructure.from_pymatgen_structure(
                diamond, scaling_matrix=scaling_matrix)
        super().__init__(structure=diamond, **kwargs)


class RocksaltStructure(Crystal3DStructure):
    """Abstract representation of caesium chloride structure."""
    def __init__(self, a=lattparams['rock_salt'], basis=['Na', 'Cl'],
                 scaling_matrix=None, **kwargs):
        rock_salt = \
            pymatgen_structure(225, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0], [0.5, 0.5, 0.5]],
                               classmethod='from_spacegroup',
                               scaling_matrix=scaling_matrix)
        rock_salt = \
            Crystal3DStructure.from_pymatgen_structure(
                rock_salt, scaling_matrix=scaling_matrix)
        super().__init__(structure=rock_salt, **kwargs)

NaClStructure = RockSaltStructure = RocksaltStructure


class ZincblendeStructure(Crystal3DStructure):
    """Abstract representation of caesium chloride structure."""
    def __init__(self, a=lattparams['zincblende'], basis=['Zn', 'Fe'],
                 scaling_matrix=None, **kwargs):
        zincblende = \
            pymatgen_structure(216, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0], [0.25, 0.25, 0.25]],
                               classmethod='from_spacegroup',
                               scaling_matrix=scaling_matrix)
        zincblende = \
            Crystal3DStructure.from_pymatgen_structure(
                zincblende, scaling_matrix=scaling_matrix)
        super().__init__(structure=zincblende, **kwargs)

SphaleriteStructure = ZincBlendeStructure = ZincblendeStructure


class CubicStructure(Crystal3DStructure):

    def __init__(self, centering, lattice=None, a=None, basis=None,
                 coords=None, scaling_matrix=None, structure=None, **kwargs):

        if lattice is None and structure is None:
            lattice = \
                Crystal3DLattice.cubic(
                    CubicStructure.get_lattice_parameter(centering, a, basis))
        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         scaling_matrix=scaling_matrix, structure=structure,
                         **kwargs)

    @classmethod
    def get_lattice_parameter(cls, centering, a, basis):
        if a is not None:
            return a

        if basis is None:
            raise ValueError('Missing required kwarg `basis`')
        elif isinstance(basis, (tuple, list)):
            if len(basis) != 1:
                raise ValueError('Expected a single element basis')
            else:
                basis = basis[0]

        if not isinstance(basis, str):
            raise ValueError('Expected `str` object for basis')

        if basis not in lattparams['cubic'][centering]:
            raise ValueError('Specify lattice constant `a` for '
                             'given basis {}'.format(basis))
        return lattparams['cubic'][centering][basis]

    @classmethod
    def from_spacegroup(cls, centering, sg, a=None, basis=None, coords=None,
                        scaling_matrix=None):
        lattice = \
            Crystal3DLattice.cubic(
                CubicStructure.get_lattice_parameter(centering, a, basis))
        if coords is None:
            coords = [[0, 0, 0]]
        return super().from_spacegroup(sg, lattice, basis, coords,
                                       scaling_matrix=scaling_matrix)


class BCCStructure(CubicStructure):
    def __init__(self, **kwargs):
        super().__init__('BCC', **kwargs)

    @classmethod
    def from_spacegroup(cls, *args, **kwargs):
        return super().from_spacegroup('BCC', *args, **kwargs)


class FCCStructure(CubicStructure):
    def __init__(self, **kwargs):
        super().__init__('FCC', **kwargs)

    @classmethod
    def from_spacegroup(cls, *args, **kwargs):
        return super().from_spacegroup('FCC', *args, **kwargs)


class Gold(FCCStructure):
    def __init__(self, **kwargs):
        kwargs['basis'] = 'Au'
        gold = FCCStructure.from_spacegroup(225, **kwargs)
        super().__init__(structure=gold, **kwargs)


class Copper(FCCStructure):
    def __init__(self, **kwargs):
        kwargs['basis'] = 'Cu'
        copper = FCCStructure.from_spacegroup(225, **kwargs)
        super().__init__(structure=copper, **kwargs)


class HexagonalClosePackedStructure(Crystal3DStructure):
    """Abstract representation of hexagonal close-packed structure."""
    pass

HCPStructure = HexagonalClosePackedStructure


class CubicClosePackedStructure(Crystal3DStructure):
    """Abstract representation of cubic close-packed structure."""
    pass

CCPStructure = CubicClosePackedStructure


class HexagonalStructure(Crystal3DStructure):
    @classmethod
    def from_spacegroup(cls, sg, a, c, basis, coords, scaling_matrix=None):
        lattice = Crystal3DLattice.hexagonal(a, c)
        return super().from_spacegroup(sg, lattice, basis, coords,
                                       scaling_matrix=scaling_matrix)


class AlphaQuartz(HexagonalStructure):
    def __init__(self, a=lattparams['alpha_quartz']['a'],
                 c=lattparams['alpha_quartz']['c'], scaling_matrix=None,
                 **kwargs):
        lattice = Crystal3DLattice.hexagonal(a, c)
        basis = Atoms(3 * ["Si"] + 6 * ["O"])
        coords = [[0.4697, 0.0000, 0.0000],
                  [0.0000, 0.4697, 0.6667],
                  [0.5305, 0.5303, 0.3333],
                  [0.4133, 0.2672, 0.1188],
                  [0.2672, 0.4133, 0.5479],
                  [0.7328, 0.1461, 0.7855],
                  [0.5867, 0.8539, 0.2145],
                  [0.8539, 0.5867, 0.4521],
                  [0.1461, 0.7328, 0.8812]]
        alpha_quartz = pymatgen_structure(lattice.cell_matrix,
                                          basis.symbols, coords,
                                          scaling_matrix=scaling_matrix)
        alpha_quartz = \
            Crystal3DStructure.from_pymatgen_structure(
                alpha_quartz, scaling_matrix=scaling_matrix)

        # alpha_quartz = \
        #     HexagonalStructure.from_spacegroup(154, a, c, ["Si", "O"],
        #                                        [[0.4697, 0.0000, 0.0000],
        #                                         [0.4135, 0.2669, 0.1191]],
        #                                        scaling_matrix=scaling_matrix)

        super().__init__(structure=alpha_quartz, **kwargs)
