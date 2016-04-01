# -*- coding: utf-8 -*-
"""
===============================================================================
Crystal structure classes (:mod:`sknano.core.structures._xtal_structures`)
===============================================================================

.. currentmodule:: sknano.core.structures._xtal_structures

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.atoms import BasisAtom, BasisAtoms
from sknano.core.crystallography import Crystal2DLattice, Crystal3DLattice
from sknano.core.math import Vector
from sknano.core.refdata import lattice_parameters
from ._base import CrystalStructureBase
from ._extras import pymatgen_structure  # , supercell_lattice_points

__all__ = ['Crystal2DStructure', 'Crystal3DStructure',
           'CaesiumChlorideStructure', 'CsClStructure',
           'DiamondStructure',
           'RocksaltStructure', 'RockSaltStructure', 'NaClStructure',
           'SphaleriteStructure', 'ZincblendeStructure', 'ZincBlendeStructure',
           'CubicStructure', 'BCCStructure', 'FCCStructure',
           'Iron', 'Copper', 'Gold',
           'CubicClosePackedStructure', 'CCPStructure',
           'HexagonalClosePackedStructure', 'HCPStructure',
           'HexagonalStructure', 'AlphaQuartz', 'BetaQuartz', 'MoS2']

cubic_lattice_parameters = lattice_parameters['cubic']


class Crystal2DStructure(CrystalStructureBase):
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
                 scaling_matrix=None, structure=None, **kwargs):

        if structure is not None:
            lattice = structure.lattice
            basis = structure.basis

        if not isinstance(lattice, Crystal2DLattice):
            lattice = Crystal2DLattice(cell_matrix=lattice)

        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         cartesian=cartesian, scaling_matrix=scaling_matrix,
                         **kwargs)

    @classmethod
    def from_pymatgen_structure(cls, structure):
        """Classmethod to generate structure using pymatgen."""
        atoms = BasisAtoms()
        for site in structure.sites:
            atoms.append(BasisAtom(element=site.specie.symbol,
                                   x=site.x, y=site.y, z=site.z))
        return cls(lattice=Crystal2DLattice(
                   cell_matrix=structure.lattice.matrix), basis=atoms)

    @classmethod
    def from_spacegroup(cls, sg, lattice=None, basis=None, coords=None,
                        **kwargs):
        """Return a `Crystal2DStructure` from a spacegroup number/symbol.

        Parameters
        ----------
        sg : :class:`~python:int` or :class:`~python:str`
        lattice : :class:`Crystal2DLattice`
        basis : :class:`~python:list` of :class:`~python:str`\ s
        coords : :class:`~python:list` of :class:`~python:float`\ s

        Returns
        -------
        :class:`Crystal2DStructure`

        Notes
        -----
        Under the hood this method first creates a pymatgen
        :class:`~pymatgen:Structure`

        See Also
        --------
        pymatgen_structure

        """
        if not isinstance(basis, list):
            basis = [basis]
        if len(basis) != len(coords):
            if isinstance(coords, list) and len(coords) != 0 and \
                    isinstance(coords[0], (int, float)):
                coords = [coords]
                cls.from_spacegroup(sg, lattice, basis, coords)

        structure = \
            pymatgen_structure(sg, lattice.cell_matrix, basis, coords,
                               classmethod='from_spacegroup')
        return cls.from_pymatgen_structure(structure)


class Crystal3DStructure(CrystalStructureBase):
    """Base class for 3D crystal structures.

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
                 scaling_matrix=None, structure=None, **kwargs):

        if structure is not None:
            lattice = structure.lattice
            basis = structure.basis

        if not isinstance(lattice, Crystal3DLattice):
            lattice = Crystal3DLattice(cell_matrix=lattice)

        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         cartesian=cartesian, scaling_matrix=scaling_matrix,
                         **kwargs)

    @classmethod
    def from_pymatgen_structure(cls, structure):
        """Return a `Crystal3DStructure` from a \
            :class:`pymatgen:pymatgen.core.Structure`.

        Parameters
        ----------
        structure : :class:`pymatgen:pymatgen.core.Structure`

        Returns
        -------
        :class:`Crystal3DStructure`

        """
        atoms = BasisAtoms()
        for site in structure.sites:
            atoms.append(BasisAtom(site.specie.symbol,
                                   x=site.x, y=site.y, z=site.z))
        return cls(lattice=Crystal3DLattice(
                   cell_matrix=structure.lattice.matrix), basis=atoms)

    @classmethod
    def from_spacegroup(cls, sg, lattice=None, basis=None, coords=None,
                        structure=None, **kwargs):
        """Return a `Crystal3DStructure` from a spacegroup number/symbol.

        Parameters
        ----------
        sg : :class:`~python:int` or :class:`~python:str`
        lattice : :class:`Crystal3DLattice`
        basis : :class:`~python:list` of :class:`~python:str`\ s
        coords : :class:`~python:list` of :class:`~python:float`\ s

        Returns
        -------
        :class:`Crystal3DStructure`

        Notes
        -----
        Under the hood, this method generates a :class:`Crystal3DStructure`
        from a :class:`~pymatgen:pymatgen.core.Structure`.

        See Also
        --------
        pymatgen_structure

        """
        if structure is not None and coords is not None:
            basis = structure.basis.copy()
            lattice = basis.lattice
            basis_elements = list(set(basis.elements))
            if len(basis_elements) != len(coords):
                raise ValueError('You must specify coordinates of structure '
                                 'basis')

            xtal_basis = BasisAtoms()
            xtal_structure = cls.from_spacegroup(sg, lattice=lattice,
                                                 basis=basis_elements,
                                                 coords=coords)
            basis.center_centroid()
            for tvec in xtal_structure.basis.r:
                basis_copy = basis.copy()
                basis_copy.rotate(axis=Vector(np.random.rand(3)).unit_vector,
                                  angle=2*np.pi*np.random.rand(1))
                basis_copy.translate(tvec)
                xtal_basis.extend(basis_copy)

            # for atom in basis:
            #     for tvec in xtal_structure.basis.r:
            #         xs, ys, zs = \
            #             lattice.cartesian_to_fractional(atom.r + tvec)
            #         xtal_basis.append(BasisAtom(atom.element,
            #                                     lattice=lattice,
            #                                     xs=xs, ys=ys, zs=zs))

            xtal_structure.basis = xtal_basis
            return xtal_structure

        if not isinstance(basis, list):
            basis = [basis]
        if len(basis) != len(coords):
            if isinstance(coords, list) and len(coords) != 0 and \
                    isinstance(coords[0], (int, float)):
                coords = [coords]
                cls.from_spacegroup(sg, lattice=lattice, basis=basis,
                                    coords=coords, **kwargs)

        structure = \
            pymatgen_structure(sg, lattice.cell_matrix, basis, coords,
                               classmethod='from_spacegroup')

        return cls.from_pymatgen_structure(structure, **kwargs)


class CubicStructure(Crystal3DStructure):
    """Base class for a cubic `Crystal3DStructure`.

    Parameters
    ----------
    centering : :class:`~python:str`
    lattice : :class:`Crystal3DLattice`, optional
    a : float, optional
    basis : :class:`~python:list`, optional
    coords : :class:`~python:list`, optional
    scaling_matrix : :class:`~python:int` or :class:`~python:list`, optional
    structure : :class:`Crystal3DStructure`, optional

    """
    def __init__(self, *args, centering=None, lattice=None, basis=None,
                 coords=None, scaling_matrix=None, structure=None,
                 **kwargs):

        if len(args) == 1 and basis is None:
            basis = args[0]

        if lattice is None and structure is None:
            lattice = CubicStructure.get_lattice(basis=basis,
                                                 centering=centering,
                                                 **kwargs)
        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         scaling_matrix=scaling_matrix, structure=structure,
                         **kwargs)

    @classmethod
    def get_lattice(cls, name=None, centering=None, basis=None, **kwargs):
        """Return cubic lattice."""
        if centering is None or not isinstance(centering, str):
            raise TypeError('Expected str for `centering` parameter')
        elif centering not in cubic_lattice_parameters:
            raise ValueError('Unrecognized `centering` value')

        a = kwargs.pop('a',
                       cubic_lattice_parameters[centering].get(name, None))
        if a is None:
            if basis is None:
                raise ValueError('`basis` cannot be None')
            elif isinstance(basis, (tuple, list)):
                if len(basis) != 1:
                    raise ValueError('Expected a single element basis')
                else:
                    basis = basis[0]

            if not isinstance(basis, str):
                raise ValueError('Expected `str` object for basis')

            if basis not in cubic_lattice_parameters[centering]:
                raise ValueError('No lattice parameters for given basis '
                                 '{}'.format(basis))

            a = cubic_lattice_parameters[centering][basis]
        lattice = Crystal3DLattice.cubic(a)
        return lattice

    @classmethod
    def from_spacegroup(cls, *args, lattice=None, basis=None, coords=None,
                        structure=None, centering=None, **kwargs):
        """Return :class:`CubicStructure` from spacegroup."""
        if len(args) == 2 and basis is None:
            sg, basis = args

        if len(args) == 1:
            sg = args[0]

        if lattice is None and structure is None:
            lattice = CubicStructure.get_lattice(basis=basis,
                                                 centering=centering,
                                                 **kwargs)
        if coords is None:
            coords = [[0, 0, 0]]
        return super().from_spacegroup(sg, lattice=lattice, basis=basis,
                                       coords=coords, structure=structure)


class CaesiumChlorideStructure(CubicStructure):
    """Abstract representation of caesium chloride structure."""
    def __init__(self, lattice=None, basis=['Cs', 'Cl'], **kwargs):

        if lattice is None:
            lattice = self.get_lattice(name='caesium_chloride',
                                       centering='cP', **kwargs)

        caesium_chloride = \
            pymatgen_structure(221, lattice.cell_matrix, basis,
                               [[0, 0, 0], [0.5, 0.5, 0.5]],
                               classmethod='from_spacegroup')
        caesium_chloride = \
            Crystal3DStructure.from_pymatgen_structure(caesium_chloride)

        super().__init__(structure=caesium_chloride, **kwargs)

CsClStructure = CaesiumChlorideStructure


class DiamondStructure(CubicStructure):
    """Abstract representation of diamond structure."""
    def __init__(self, lattice=None, basis=['C'], **kwargs):

        if lattice is None:
            lattice = self.get_lattice(name='diamond', centering='FCC',
                                       **kwargs)

        diamond = \
            pymatgen_structure(227, lattice.cell_matrix, basis, [[0, 0, 0]],
                               classmethod='from_spacegroup')
        diamond = \
            Crystal3DStructure.from_pymatgen_structure(diamond)
        super().__init__(structure=diamond, **kwargs)


class RocksaltStructure(CubicStructure):
    """Abstract representation of rock salt structure."""
    def __init__(self, lattice=None, basis=['Na', 'Cl'], **kwargs):
        if lattice is None:
            lattice = self.get_lattice(name='rocksalt', centering='FCC',
                                       **kwargs)

        rock_salt = \
            pymatgen_structure(225, lattice.cell_matrix, basis,
                               [[0, 0, 0], [0.5, 0.5, 0.5]],
                               classmethod='from_spacegroup')
        rock_salt = \
            Crystal3DStructure.from_pymatgen_structure(rock_salt)
        super().__init__(structure=rock_salt, **kwargs)

NaClStructure = RockSaltStructure = RocksaltStructure


class ZincblendeStructure(CubicStructure):
    """Abstract representation of zinc blende structure."""
    def __init__(self, lattice=None, basis=['Zn', 'S'], **kwargs):
        if lattice is None:
            lattice = self.get_lattice(name='zincblende', centering='FCC',
                                       **kwargs)

        zincblende = \
            pymatgen_structure(216, lattice.cell_matrix, basis,
                               [[0, 0, 0], [0.25, 0.25, 0.25]],
                               classmethod='from_spacegroup')
        zincblende = \
            Crystal3DStructure.from_pymatgen_structure(zincblende)
        super().__init__(structure=zincblende, **kwargs)

SphaleriteStructure = ZincBlendeStructure = ZincblendeStructure


class BCCStructure(CubicStructure):
    """BCC structure class."""
    def __init__(self, *args, **kwargs):
        kwargs['centering'] = 'bcc'
        kwargs['structure'] = \
            CubicStructure.from_spacegroup(229, *args, **kwargs)
        super().__init__(*args, **kwargs)


class FCCStructure(CubicStructure):
    """FCC structure class."""
    def __init__(self, *args, **kwargs):
        kwargs['centering'] = 'fcc'
        kwargs['structure'] = \
            CubicStructure.from_spacegroup(225, *args, **kwargs)
        super().__init__(*args, **kwargs)


class Iron(BCCStructure):
    """Iron structure."""
    def __init__(self, **kwargs):
        kwargs['basis'] = 'Fe'
        super().__init__(**kwargs)


class Gold(FCCStructure):
    """Gold structure."""
    def __init__(self, **kwargs):
        kwargs['basis'] = 'Au'
        super().__init__(**kwargs)


class Copper(FCCStructure):
    """Copper structure."""
    def __init__(self, **kwargs):
        kwargs['basis'] = 'Cu'
        super().__init__(**kwargs)


class HexagonalClosePackedStructure(Crystal3DStructure):
    """Abstract representation of hexagonal close-packed structure."""
    pass

HCPStructure = HexagonalClosePackedStructure


class CubicClosePackedStructure(Crystal3DStructure):
    """Abstract representation of cubic close-packed structure."""
    pass

CCPStructure = CubicClosePackedStructure


class HexagonalStructure(Crystal3DStructure):
    """Base :class:`Crystal3DStructure` for a hexagonal crystal system.

    Parameters
    ----------
    centering : :class:`~python:str`
    lattice : :class:`Crystal3DLattice`, optional
    a, c : :class:`~python:float`, optional
    basis : :class:`~python:list`, optional
    coords : :class:`~python:list`, optional
    scaling_matrix : :class:`~python:int` or :class:`~python:list`, optional
    structure : :class:`Crystal3DStructure`, optional

    """
    def __init__(self, *args, lattice=None, basis=None, coords=None,
                 scaling_matrix=None, structure=None, **kwargs):

        if len(args) == 1 and basis is None:
            basis = args[0]

        if lattice is None and structure is None:
            lattice = self.get_lattice(basis=basis, **kwargs)

        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         scaling_matrix=scaling_matrix, structure=structure,
                         **kwargs)

    @classmethod
    def get_lattice(cls, name=None, basis=None, **kwargs):
        """Return hexagonal :class:`Crystal3DLattice`."""
        a, c = [kwargs.pop(k, default) for k, default in
                zip(('a', 'c'),
                    lattice_parameters['hexagonal'].get(name, (None, None)))]

        if a is None or c is None:
            if basis is None:
                raise ValueError('`basis` cannot be None.')
            elif isinstance(basis, (tuple, list)):
                if len(basis) != 1:
                    raise ValueError('Expected a single element basis')
                else:
                    basis = basis[0]

            if not isinstance(basis, str):
                raise ValueError('Expected `str` object for basis')

            if basis not in lattice_parameters['hexagonal']:
                raise ValueError('No lattice constants for '
                                 'given basis {}'.format(basis))

            a, c = lattice_parameters['hexagonal'][basis]

        lattice = Crystal3DLattice.hexagonal(a, c)
        return lattice

    @classmethod
    def from_spacegroup(cls, *args, lattice=None, basis=None, coords=None,
                        structure=None, **kwargs):
        """Return :class:`HexagonalStructure` from space group."""
        if len(args) == 2 and basis is None:
            sg, basis = args

        if len(args) == 1:
            sg = args[0]

        if lattice is None and structure is None:
            lattice = HexagonalStructure.get_lattice(basis=basis, **kwargs)

        if coords is None:
            coords = [[0, 0, 0]]
        return super().from_spacegroup(sg, lattice=lattice, basis=basis,
                                       coords=coords, structure=structure)


class MoS2(HexagonalStructure):
    """Molybdenum disulphide structure class."""
    def __init__(self, lattice=None, **kwargs):

        if lattice is None:
            lattice = self.get_lattice(name='MoS2', **kwargs)
        basis = ['Mo', 'S']
        molybdenum_disulphide = \
            pymatgen_structure(194, lattice.cell_matrix, basis,
                               [[1/3, 2/3, 1/4], [1/3, 2/3, 0.621]],
                               classmethod='from_spacegroup')
        molybdenum_disulphide = \
            Crystal3DStructure.from_pymatgen_structure(molybdenum_disulphide)

        super().__init__(structure=molybdenum_disulphide, **kwargs)


class AlphaQuartz(HexagonalStructure):
    """Alpha quartz structure class."""
    def __init__(self, lattice=None, **kwargs):
        if lattice is None:
            lattice = self.get_lattice(name='alpha_quartz', **kwargs)
        # basis = BasisAtoms(3 * ["Si"] + 6 * ["O"])
        basis = 3 * ['Si'] + 6 * ['O']
        coords = [[0.4697, 0.0000, 0.0000],
                  [0.0000, 0.4697, 0.6667],
                  [0.5305, 0.5303, 0.3333],
                  [0.4133, 0.2672, 0.1188],
                  [0.2672, 0.4133, 0.5479],
                  [0.7328, 0.1461, 0.7855],
                  [0.5867, 0.8539, 0.2145],
                  [0.8539, 0.5867, 0.4521],
                  [0.1461, 0.7328, 0.8812]]
        alpha_quartz = pymatgen_structure(lattice.cell_matrix, basis, coords)
        alpha_quartz = \
            Crystal3DStructure.from_pymatgen_structure(alpha_quartz)

        # alpha_quartz = \
        #     HexagonalStructure.from_spacegroup(154, lattice, ["Si", "O"],
        #                                        [[0.4697, 0.0000, 0.0000],
        #                                         [0.4135, 0.2669, 0.1191]])

        super().__init__(structure=alpha_quartz, **kwargs)


class BetaQuartz(HexagonalStructure):
    """Beta quartz structure class."""
    def __init__(self, lattice=None, **kwargs):
        if lattice is None:
            lattice = self.get_lattice(name='beta_quartz', **kwargs)

        basis = ["Si", "O"]
        # basis = 3 * ['Si'] + 6 * ['O']

        # basis = BasisAtoms(3 * ["Si"] + 6 * ["O"])
        # coords = [[0.4697, 0.0000, 0.0000],
        #           [0.0000, 0.4697, 0.6667],
        #           [0.5305, 0.5303, 0.3333],
        #           [0.4133, 0.2672, 0.1188],
        #           [0.2672, 0.4133, 0.5479],
        #           [0.7328, 0.1461, 0.7855],
        #           [0.5867, 0.8539, 0.2145],
        #           [0.8539, 0.5867, 0.4521],
        #           [0.1461, 0.7328, 0.8812]]
        # beta_quartz = pymatgen_structure(lattice.cell_matrix,
        #                                  basis.symbols, coords)
        # alpha_quartz = \
        #     Crystal3DStructure.from_pymatgen_structure(alpha_quartz)

        beta_quartz = pymatgen_structure(180, lattice.cell_matrix, basis,
                                         [[0.5000, 0.0000, 0.0000],
                                          [0.4147, 0.2078, 0.1667]],
                                         classmethod='from_spacegroup')
        beta_quartz = Crystal3DStructure.from_pymatgen_structure(beta_quartz)

        super().__init__(structure=beta_quartz, **kwargs)
