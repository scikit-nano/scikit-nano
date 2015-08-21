# -*- coding: utf-8 -*-
"""
===============================================================================
Crystal structure classes (:mod:`sknano.core.crystallography._xtal_structures`)
===============================================================================

.. currentmodule:: sknano.core.crystallography._xtal_structures

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering

# import numpy as np

from sknano.core import BaseClass
from sknano.core.atoms import BasisAtom as Atom, BasisAtoms as Atoms, \
    StructureAtoms
from sknano.core.refdata import lattice_parameters as lattparams

from ._extras import pymatgen_structure
from ._xtal_cells import CrystalCell, UnitCell, SuperCell
from ._xtal_lattices import Crystal2DLattice, Crystal3DLattice

__all__ = ['BaseStructureMixin', 'BaseStructure', 'StructureData',
           'CrystalStructureBase', 'Crystal2DStructure',
           'CrystalStructure', 'Crystal3DStructure',
           'CaesiumChlorideStructure', 'CsClStructure',
           'DiamondStructure',
           'RocksaltStructure', 'RockSaltStructure', 'NaClStructure',
           'SphaleriteStructure', 'ZincblendeStructure', 'ZincBlendeStructure',
           'BCCStructure', 'FCCStructure', 'Iron', 'Copper', 'Gold',
           'CubicClosePackedStructure', 'CCPStructure',
           'HexagonalClosePackedStructure', 'HCPStructure',
           'HexagonalStructure', 'AlphaQuartz', 'MoS2']


class BaseStructureMixin:
    """Mixin class for crystal structures."""

    # def __deepcopy__(self, memo):
    #     from copy import deepcopy
    #     cp = self.__class__()
    #     memo[id(self)] = cp
    #     for attr in dir(self):
    #         if not attr.startswith('_'):
    #             setattr(cp, attr, deepcopy(getattr(self, attr), memo))
    #     return cp

    def __getattr__(self, name):
        try:
            return getattr(self.atoms, name)
        except AttributeError:
            try:
                return getattr(self.crystal_cell, name)
            except AttributeError:
                return super().__getattr__(name)

    @property
    def atoms(self):
        return self._atoms

    @property
    def crystal_cell(self):
        return self._crystal_cell

    @crystal_cell.setter
    def crystal_cell(self, value):
        self._crystal_cell = value

    @property
    def basis(self):
        return self.crystal_cell.basis

    @basis.setter
    def basis(self, value):
        self.crystal_cell.basis = value

    @property
    def lattice(self):
        return self.crystal_cell.lattice

    @lattice.setter
    def lattice(self, value):
        self.crystal_cell.lattice = value

    @property
    def scaling_matrix(self):
        return self.crystal_cell.scaling_matrix

    @scaling_matrix.setter
    def scaling_matrix(self, value):
        self.crystal_cell.scaling_matrix = value

    @property
    def unit_cell(self):
        return self.crystal_cell.unit_cell

    @unit_cell.setter
    def unit_cell(self, value):
        self.crystal_cell.unit_cell = value

    @property
    def structure_data(self):
        return self

    @property
    def structure(self):
        return self

    def clear(self):
        self.atoms.clear()

    def rotate(self, **kwargs):
        """Rotate crystal cell lattice, basis, and unit cell."""
        self.crystal_cell.rotate(**kwargs)

    def translate(self, t, fix_anchor_points=True):
        """Translate crystal cell basis."""
        self.crystal_cell.translate(t, fix_anchor_points=fix_anchor_points)


class BaseStructure(BaseStructureMixin):
    """Base structure class for structure data."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._atoms = StructureAtoms()
        self._crystal_cell = CrystalCell()

StructureData = BaseStructure


@total_ordering
class CrystalStructureBase(BaseStructure, BaseClass):
    """Base class for abstract representions of crystal structures.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional
    scaling_matrix : {:class:`~python:int`, :class:`~python:list`}, optional

    """

    def __init__(self, lattice=None, basis=None, coords=None,
                 cartesian=False, scaling_matrix=None, **kwargs):

        super().__init__(**kwargs)
        self.unit_cell = UnitCell(lattice=lattice, basis=basis,
                                  coords=coords, cartesian=cartesian)

        self.scaling_matrix = scaling_matrix
        self.fmtstr = self.unit_cell.fmtstr + \
            ", scaling_matrix={scaling_matrix!r}"

    def __dir__(self):
        return dir(self.crystal_cell)

    def __eq__(self, other):
        if isinstance(other, CrystalStructureBase):
            return self is other or self.crystal_cell == other.crystal_cell

    def __lt__(self, other):
        if isinstance(other, CrystalStructureBase):
            return self.crystal_cell < other.crystal_cell

    def make_supercell(self, scaling_matrix):
        return SuperCell(self.unit_cell, scaling_matrix)

    def todict(self):
        attrdict = self.unit_cell.todict()
        attrdict.update(dict(scaling_matrix=self.scaling_matrix.tolist()))
        return attrdict


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
        atoms = Atoms()
        for site in structure.sites:
            atoms.append(Atom(element=site.specie.symbol,
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
        Under the hood this method first s a pymatgen
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
        atoms = Atoms()
        for site in structure.sites:
            atoms.append(Atom(site.specie.symbol,
                              x=site.x, y=site.y, z=site.z))
        return cls(lattice=Crystal3DLattice(
                   cell_matrix=structure.lattice.matrix), basis=atoms)

    @classmethod
    def from_spacegroup(cls, sg, lattice=None, basis=None, coords=None,
                        **kwargs):
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
        Under the hood this method generates a
        :class:`~pymatgen:pymatgen.core.Structure`

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
                cls.from_spacegroup(sg, lattice=lattice, basis=basis,
                                    coords=coords, **kwargs)

        structure = \
            pymatgen_structure(sg, lattice.cell_matrix, basis, coords,
                               classmethod='from_spacegroup')

        return cls.from_pymatgen_structure(structure, **kwargs)

CrystalStructure = Crystal3DStructure


class MoS2(Crystal3DStructure):
    """Molybdenum disulphide structure class."""
    def __init__(self, a=lattparams['molybdenum_disulphide']['a'],
                 c=lattparams['molybdenum_disulphide']['c'],
                 basis=['Mo', 'S'], **kwargs):

        molybdenum_disulphide = \
            pymatgen_structure(194,
                               Crystal3DLattice.hexagonal(a, c).cell_matrix,
                               basis, [[1/3, 2/3, 1/4], [1/3, 2/3, 0.621]],
                               classmethod='from_spacegroup')
        molybdenum_disulphide = \
            Crystal3DStructure.from_pymatgen_structure(molybdenum_disulphide)

        super().__init__(structure=molybdenum_disulphide, **kwargs)


class CaesiumChlorideStructure(Crystal3DStructure):
    """Abstract representation of caesium chloride structure."""
    def __init__(self, a=lattparams['caesium_chloride'], basis=['Cs', 'Cl'],
                 **kwargs):
        caesium_chloride = \
            pymatgen_structure(221, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0], [0.5, 0.5, 0.5]],
                               classmethod='from_spacegroup')
        caesium_chloride = \
            Crystal3DStructure.from_pymatgen_structure(caesium_chloride)

        super().__init__(structure=caesium_chloride, **kwargs)

CsClStructure = CaesiumChlorideStructure


class DiamondStructure(Crystal3DStructure):
    """Abstract representation of diamond structure."""
    def __init__(self, a=lattparams['diamond'], basis=['C'], **kwargs):
        diamond = \
            pymatgen_structure(227, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0]],
                               classmethod='from_spacegroup')
        diamond = \
            Crystal3DStructure.from_pymatgen_structure(diamond)
        super().__init__(structure=diamond, **kwargs)


class RocksaltStructure(Crystal3DStructure):
    """Abstract representation of caesium chloride structure."""
    def __init__(self, a=lattparams['rock_salt'], basis=['Na', 'Cl'],
                 **kwargs):
        rock_salt = \
            pymatgen_structure(225, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0], [0.5, 0.5, 0.5]],
                               classmethod='from_spacegroup')
        rock_salt = \
            Crystal3DStructure.from_pymatgen_structure(rock_salt)
        super().__init__(structure=rock_salt, **kwargs)

NaClStructure = RockSaltStructure = RocksaltStructure


class ZincblendeStructure(Crystal3DStructure):
    """Abstract representation of caesium chloride structure."""
    def __init__(self, a=lattparams['zincblende'], basis=['Zn', 'Fe'],
                 **kwargs):
        zincblende = \
            pymatgen_structure(216, Crystal3DLattice.cubic(a).cell_matrix,
                               basis, [[0, 0, 0], [0.25, 0.25, 0.25]],
                               classmethod='from_spacegroup')
        zincblende = \
            Crystal3DStructure.from_pymatgen_structure(zincblende)
        super().__init__(structure=zincblende, **kwargs)

SphaleriteStructure = ZincBlendeStructure = ZincblendeStructure


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
    def __init__(self, *args, a=None, centering=None, lattice=None,
                 basis=None, coords=None, scaling_matrix=None, structure=None,
                 **kwargs):

        if len(args) == 1 and basis is None:
            basis = args[0]

        self.centering = centering

        if lattice is None and structure is None:
            lattice = \
                Crystal3DLattice.cubic(
                    CubicStructure.get_lattice_parameter(
                        a=a, basis=basis, centering=centering))
        super().__init__(lattice=lattice, basis=basis, coords=coords,
                         scaling_matrix=scaling_matrix, structure=structure,
                         **kwargs)

    @classmethod
    def get_lattice_parameter(cls, a=None, basis=None, centering=None):
        if a is not None:
            return a

        if basis is None or centering is None:
            raise ValueError('\nBoth the `basis` and `centering` kwargs '
                             'are required')
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
    def from_spacegroup(cls, *args, lattice=None, basis=None, coords=None,
                        a=None, centering=None, **kwargs):
        if len(args) == 2 and basis is None:
            sg, basis = args

        if len(args) == 1:
            sg = args[0]

        if lattice is None:
            lattice = \
                Crystal3DLattice.cubic(
                    CubicStructure.get_lattice_parameter(a=a, basis=basis,
                                                         centering=centering))
        if coords is None:
            coords = [[0, 0, 0]]
        return super().from_spacegroup(sg, lattice=lattice, basis=basis,
                                       coords=coords)


class BCCStructure(CubicStructure):
    """BCC structure class."""
    def __init__(self, *args, **kwargs):
        kwargs['centering'] = 'BCC'
        structure = CubicStructure.from_spacegroup(229, *args, **kwargs)
        super().__init__(*args, structure=structure, **kwargs)


class FCCStructure(CubicStructure):
    """FCC structure class."""
    def __init__(self, *args, **kwargs):
        kwargs['centering'] = 'FCC'
        structure = CubicStructure.from_spacegroup(225, *args, **kwargs)
        super().__init__(*args, structure=structure, **kwargs)


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
    @classmethod
    def from_spacegroup(cls, sg, a, c, basis, coords):
        lattice = Crystal3DLattice.hexagonal(a, c)
        return super().from_spacegroup(sg, lattice=lattice, basis=basis,
                                       coords=coords)


class AlphaQuartz(HexagonalStructure):
    """Alpha quartz structure class."""
    def __init__(self, a=lattparams['alpha_quartz']['a'],
                 c=lattparams['alpha_quartz']['c'], **kwargs):
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
                                          basis.symbols, coords)
        alpha_quartz = \
            Crystal3DStructure.from_pymatgen_structure(alpha_quartz)

        # alpha_quartz = \
        #     HexagonalStructure.from_spacegroup(154, a, c, ["Si", "O"],
        #                                        [[0.4697, 0.0000, 0.0000],
        #                                         [0.4135, 0.2669, 0.1191]],
        #                                        scaling_matrix=scaling_matrix)

        super().__init__(structure=alpha_quartz, **kwargs)
