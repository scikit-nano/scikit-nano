# -*- coding: utf-8 -*-
"""
==============================================================================
Base structure classes (:mod:`sknano.core.structures.base`)
==============================================================================

.. currentmodule:: sknano.core.structures.base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering

import numpy as np

from sknano.core import BaseClass, TabulateMixin
from sknano.core.atoms import Atoms, BasisAtoms, StructureAtoms, \
    vdw_radius_from_basis
from sknano.core.crystallography import CrystalCell, UnitCell, SuperCell
from sknano.core.math import Vector
from sknano.core.refdata import aCC, element_data

from .selections import StructureSelectionMixin


__all__ = ['StructureMixin', 'StructureBase',
           'CrystalStructureBase', 'NanoStructureBase']

r_CC_vdw = element_data['C']['VanDerWaalsRadius']

_list_methods = ('append', 'extend', 'insert', 'remove', 'pop', 'clear',
                 'index', 'count', 'sort', 'reverse', 'copy')


class StructureMixin:
    """Mixin class for crystal structures."""
    def __getattr__(self, name):
        try:
            if name not in _list_methods and not name.startswith('_') \
                    and self.atoms is None:
                raise ValueError
            return getattr(self.atoms, name)
        except (AttributeError, ValueError):
            try:
                return getattr(self.crystal_cell, name)
            except AttributeError:
                return super().__getattr__(name)

    @property
    def atoms(self):
        """Structure :class:`~sknano.core.atoms.Atoms`."""
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        if not isinstance(value, Atoms):
            errmsg = '{} is not a valid `Atoms` object.'.format(value)
            raise ValueError(errmsg)
        self._atoms = value

    @property
    def crystal_cell(self):
        """Structure :class:`~sknano.core.crystallography.CrystalCell`."""
        return self._crystal_cell

    @crystal_cell.setter
    def crystal_cell(self, value):
        self._crystal_cell = value

    @property
    def basis(self):
        """Structure :class:`~sknano.core.atoms.BasisAtoms`."""
        return self.crystal_cell.basis

    @basis.setter
    def basis(self, value):
        if not isinstance(value, BasisAtoms):
            raise ValueError('Expected a `BasisAtoms` object')
        self.crystal_cell.basis = value

    @property
    def lattice(self):
        """Structure :class:`~sknano.core.crystallography.Crystal3DLattice`."""
        return self.crystal_cell.lattice

    @lattice.setter
    def lattice(self, value):
        self.crystal_cell.lattice = value
        self.atoms.lattice = self.crystal_cell.lattice

    @property
    def lattice_shift(self):
        """Lattice displacement vector."""
        return self._lattice_shift

    @lattice_shift.setter
    def lattice_shift(self, value):
        self._lattice_shift = Vector(value)

    @property
    def scaling_matrix(self):
        """:attr:`CrystalCell.scaling_matrix`."""
        return self.crystal_cell.scaling_matrix

    @scaling_matrix.setter
    def scaling_matrix(self, value):
        self.crystal_cell.scaling_matrix = value

    @property
    def unit_cell(self):
        """Structure :class:`~sknano.core.crystallography.UnitCell`."""
        return self.crystal_cell.unit_cell

    @unit_cell.setter
    def unit_cell(self, value):
        self.crystal_cell.unit_cell = value

    @property
    def structure(self):
        """An alias to `self`."""
        return self

    def clear(self):
        """Clear list of :attr:`StructureMixin.atoms`."""
        self.atoms.clear()

    def make_supercell(self, scaling_matrix, wrap_coords=False):
        """Make supercell."""
        return SuperCell(self.unit_cell, scaling_matrix,
                         wrap_coords=wrap_coords)

    def rotate(self, **kwargs):
        """Rotate crystal cell lattice, basis, and unit cell."""
        self.crystal_cell.rotate(**kwargs)
        self.atoms.rotate(**kwargs)
        # self.lattice = self.crystal_cell.lattice

    def translate(self, t, fix_anchor_points=True):
        """Translate crystal cell lattice, basis, and unit cell."""
        self.crystal_cell.translate(t, fix_anchor_points=fix_anchor_points)
        self.atoms.translate(t, fix_anchor_points=fix_anchor_points)

    def transform_lattice(self, scaling_matrix, wrap_coords=False, pbc=None):
        """Transform structure lattice."""
        if self.lattice is None:
            return
        self.lattice = self.lattice.__class__(
            cell_matrix=np.asmatrix(scaling_matrix) * self.lattice.matrix)

        if wrap_coords:
            self.crystal_cell.basis.wrap_coords(pbc=pbc)
            self.atoms.wrap_coords(pbc=pbc)

        # tvecs = \
        #     np.asarray(
        #         np.asmatrix(supercell_lattice_points(scaling_matrix)) *
        #         self.lattice.matrix)

        # if self.crystal_cell.basis is not None:
        #     basis = self.crystal_cell.basis[:]
        #     self.crystal_cell.basis = BasisAtoms()
        #     for atom in basis:
        #         xs, ys, zs = \
        #             self.lattice.cartesian_to_fractional(atom.r)
        #         if wrap_coords:
        #             xs, ys, zs = \
        #                 self.lattice.wrap_fractional_coordinate([xs, ys, zs])
        #         self.crystal_cell.basis.append(
        #             BasisAtom(atom.element, lattice=self.lattice,
        #                       xs=xs, ys=ys, zs=zs))

    # def __deepcopy__(self, memo):
    #     from copy import deepcopy
    #     cp = self.__class__()
    #     memo[id(self)] = cp
    #     for attr in dir(self):
    #         if not attr.startswith('_'):
    #             setattr(cp, attr, deepcopy(getattr(self, attr), memo))
    #     return cp


class StructureBase(StructureSelectionMixin, TabulateMixin, StructureMixin):
    """Base structure class for structure data."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._atoms = StructureAtoms()
        self._crystal_cell = CrystalCell()
        self._lattice_shift = Vector()
        self._region = None

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        atoms = self.atoms
        xtal_cell = self.crystal_cell

        title = '.'.join((objstr, atoms.__class__.__name__))
        strrep = '\n'.join((strrep, title, str(atoms)))

        title = '.'.join((objstr, xtal_cell.__class__.__name__))
        strrep = '\n'.join((strrep, title, str(xtal_cell)))
        return strrep


@total_ordering
class CrystalStructureBase(StructureBase, BaseClass):
    """Base class for abstract representions of crystal structures.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional
    scaling_matrix : {:class:`~python:int`, :class:`~python:list`}, optional

    """
    def __init__(self, *args, lattice=None, basis=None, coords=None,
                 cartesian=False, scaling_matrix=None, **kwargs):

        super().__init__(*args, **kwargs)
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

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attrdict = self.unit_cell.todict()
        attrdict.update(dict(scaling_matrix=self.scaling_matrix.tolist()))
        return attrdict


class NanoStructureBase(StructureBase, BaseClass):
    """Base class for creating abstract representations of nanostructure.

    Parameters
    ----------
    basis : {:class:`~python:str`, :class:`~python:int`, list}, optional
        Element chemical symbols or atomic numbers of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2

    """
    def __init__(self, *args, basis=None, bond=None, **kwargs):

        if basis is None:
            basis = ['C', 'C']

        if isinstance(basis, list):
            basis = basis[:]

            [basis.__setitem__(i, kwargs.pop(element)) for i, element in
             enumerate(('element1', 'element2')) if element in kwargs]

        vdw_radius = None
        if 'vdw_spacing' in kwargs:
            vdw_radius = kwargs.pop('vdw_spacing') / 2
        elif 'vdw_radius' in kwargs:
            vdw_radius = kwargs.pop('vdw_radius')

        if vdw_radius is None:
            vdw_radius = r_CC_vdw

        if bond is None:
            bond = aCC

        super().__init__(*args, **kwargs)

        self.bond = bond
        self.basis = basis
        self.vdw_radius = vdw_radius
        self.fmtstr = "bond={bond!r}, basis={basis!r}, " + \
            "vdw_radius={vdw_radius!r}"

    @property
    def basis(self):
        """:class:`NanoStructureBase` basis objects."""
        return self._basis

    @basis.setter
    def basis(self, value):
        self._basis = value
        try:
            [self.crystal_cell.update_basis(obj, index=i, step=2) for
             i, obj in enumerate(self.basis)]
        except AttributeError:
            pass

    @basis.deleter
    def basis(self):
        del self._basis

    @property
    def vdw_radius(self):
        """Van der Waals radius"""
        if self._vdw_radius is not None:
            return self._vdw_radius
        else:
            return vdw_radius_from_basis(self.basis[0], self.basis[1])

    @vdw_radius.setter
    def vdw_radius(self, value):
        self._vdw_radius = value

    @property
    def vdw_distance(self):
        """Van der Waals distance."""
        return 2 * self.vdw_radius

    @property
    def element1(self):
        """Basis element 1"""
        return self.basis[0]

    @element1.setter
    def element1(self, value):
        self.basis[0] = value
        try:
            self.crystal_cell.update_basis(value, index=0, step=2)
        except AttributeError:
            pass

    @property
    def element2(self):
        """Basis element 2"""
        return self.basis[1]

    @element2.setter
    def element2(self, value):
        self.basis[1] = value
        try:
            self.crystal_cell.update_basis(value, index=1, step=2)
        except AttributeError:
            pass

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(bond=self.bond, basis=self.basis,
                    vdw_radius=self.vdw_radius)
