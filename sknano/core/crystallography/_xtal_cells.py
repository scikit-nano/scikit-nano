# -*- coding: utf-8 -*-
"""
===========================================================================
Crystal cell classes (:mod:`sknano.core.crystallography._xtal_cells`)
===========================================================================

.. currentmodule:: sknano.core.crystallography._xtal_cells

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering
# import copy
import numbers

import numpy as np

from sknano.core import BaseClass
from sknano.core.atoms import BasisAtom, BasisAtoms
from ._extras import supercell_lattice_points

__all__ = ['CrystalCell', 'UnitCell', 'SuperCell']


@total_ordering
class UnitCell(BaseClass):
    """Base class for abstract representations of crystallographic unit cells.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional
    wrap_coords : {:class:`~python:bool`}, optional

    """

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False,
                 wrap_coords=False):

        super().__init__()

        if basis is None:
            basis = BasisAtoms()
        else:
            basis = BasisAtoms(basis)
            basis.lattice = lattice
            if lattice is not None and coords is not None:
                for atom, pos in zip(basis, coords):
                    atom.lattice = lattice
                    if not cartesian:
                        atom.rs = pos
                    else:
                        atom.rs = lattice.cartesian_to_fractional(pos)

        self.lattice = lattice
        self.basis = basis
        self.wrap_coords = wrap_coords
        self.fmtstr = "{lattice!r}, {basis!r}, {coords!r}, " + \
            "cartesian=False, wrap_coords={wrap_coords!r}"

    def __dir__(self):
        return ['lattice', 'basis']

    def __eq__(self, other):
        return self is other or \
            all([(getattr(self, attr) == getattr(other, attr)) for attr
                 in dir(self)])

    def __lt__(self, other):
        return (self.lattice < other.lattice and self.basis <= other.basis) \
            or (self.lattice <= other.lattice and self.basis < other.basis)

    def __getattr__(self, name):
        try:
            return getattr(self.lattice, name)
        except AttributeError:
            try:
                return getattr(self.basis, name)
            except AttributeError:
                return super().__getattr__(name)

    def __iter__(self):
        return iter(self.basis)

    # @property
    # def basis(self):
    #     return self._basis

    # @basis.setter
    # def basis(self, value):
    #     lattice = self.lattice
    #     coords = self.coords
    #     if value is None:
    #         value = BasisAtoms()
    #     elif value is not None:
    #         value = BasisAtoms(value, lattice=lattice)
    #         if coords is not None:
    #             for atom, pos in zip(basis, coords):
    #                 atom.lattice = lattice
    #                 if not cartesian:
    #                     atom.rs = pos
    #                 else:
    #                     atom.rs = lattice.cartesian_to_fractional(pos)

    def rotate(self, **kwargs):
        """Rotate unit cell lattice vectors and basis."""
        self.lattice.rotate(**kwargs)
        self.basis.rotate(**kwargs)

    def translate(self, t, fix_anchor_points=True):
        """Translate unit cell basis."""
        if not fix_anchor_points:
            self.lattice.translate(t)
        self.basis.translate(t, fix_anchor_points=fix_anchor_points)

    def todict(self):
        """Return `dict` of `UnitCell` parameters."""
        return dict(lattice=self.lattice, basis=self.basis.symbols.tolist(),
                    coords=self.basis.rs.tolist(),
                    wrap_coords=self.wrap_coords)


@total_ordering
class CrystalCell(BaseClass):
    """Class representation of crystal structure cell.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional
    wrap_coords : {:class:`~python:bool`}, optional
    unit_cell : :class:`~sknano.core.crystallography.UnitCell`
    scaling_matrix : {:class:`~python:int`, :class:`~python:list`}

    """

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False,
                 wrap_coords=False, unit_cell=None, scaling_matrix=None):
        super().__init__()

        if unit_cell is None and basis is not None:
            basis = BasisAtoms(basis)
            basis.lattice = lattice
            if lattice is not None and coords is not None:
                for atom, pos in zip(basis, coords):
                    atom.lattice = lattice
                    if not cartesian:
                        atom.rs = pos
                    else:
                        atom.rs = lattice.cartesian_to_fractional(pos)

        # if basis is None:
        #     basis = BasisAtoms()

        # These attributes may be reset in the `@scaling_matrix.setter`
        # method and so they need to be initialized *before* setting
        # `self.scaling_matrix`.
        self.basis = basis
        self.lattice = lattice

        self.unit_cell = unit_cell
        self.wrap_coords = wrap_coords
        self.scaling_matrix = scaling_matrix

        self.fmtstr = \
            "lattice={lattice!r}, basis={basis!r}, coords={coords!r}, " + \
            "cartesian=False, wrap_coords={wrap_coords!r}, " + \
            "unit_cell={unit_cell!r}, scaling_matrix={scaling_matrix!r}"

    def __dir__(self):
        return ['lattice', 'basis', 'unit_cell', 'scaling_matrix']

    def __eq__(self, other):
        if all([attr is not None for attr in
                (self.scaling_matrix, self.unit_cell,
                 other.scaling_matrix, other.unit_cell)]) and \
                self.scaling_matrix.shape == other.scaling_matrix.shape:
            return self is other or \
                (self.unit_cell == other.unit_cell and
                 np.allclose(self.scaling_matrix, other.scaling_matrix))
        elif all([cell is not None for cell in
                  (self.unit_cell, other.unit_cell)]):
            return self is other or \
                (self.unit_cell == other.unit_cell and
                 all([mat is None for mat in
                      (self.scaling_matrix, other.scaling_matrix)]))

    def __lt__(self, other):
        return (self.unit_cell < other.unit_cell and
                self.scaling_matrix <= other.scaling_matrix) \
            or (self.unit_cell <= other.unit_cell and
                self.scaling_matrix < other.scaling_matrix)

    def __iter__(self):
        return iter(self.basis)

    def __getattr__(self, name):
        if name != 'lattice' and self.lattice is not None:
            try:
                return getattr(self.lattice, name)
            except AttributeError:
                pass
        if name != 'basis' and self.basis.Natoms != 0:
            try:
                return getattr(self.basis, name)
            except AttributeError:
                pass
        try:
            return getattr(self.unit_cell, name)
        except AttributeError:
            return super().__getattr__(name)

    @property
    def basis(self):
        return self._basis

    @basis.setter
    def basis(self, value):
        self._basis = value
        # if self.unit_cell is not None:
        #     self.unit_cell.basis[:] = \
        #         self.basis[:self.unit_cell.basis.Natoms]

    @property
    def lattice(self):
        return self._lattice

    @lattice.setter
    def lattice(self, value):
        self._lattice = value
        if self.basis is not None:
            self.basis.lattice = self.lattice

    @property
    def unit_cell(self):
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, value):
        if value is not None and not isinstance(value, UnitCell):
            raise ValueError('Expected a `UnitCell` object')
        self._unit_cell = value
        if value is not None:
            if self.lattice is None:
                self._lattice = self.unit_cell.lattice
            if self.basis is None or self.basis.Natoms == 0:
                self._basis = self.unit_cell.basis

    @property
    def scaling_matrix(self):
        """Scaling matrix."""
        return self._scaling_matrix

    @scaling_matrix.setter
    def scaling_matrix(self, value):
        if value is None:
            self._scaling_matrix = np.asmatrix(np.ones(3, dtype=int))
            return

        if not isinstance(value, (int, float, tuple, list, np.ndarray)):
            return

        if isinstance(value, np.ndarray) and \
                ((value.shape == np.ones(3).shape and
                  np.allclose(value, np.ones(3))) or
                 (value.shape == np.eye(3).shape and
                  np.allclose(value, np.eye(3)))):
            self._scaling_matrix = np.asmatrix(value)
            return

        if isinstance(value, numbers.Number):
            value = self.lattice.nd * [int(value)]

        scaling_matrix = np.asmatrix(value, dtype=int)
        # scaling_matrix = np.asmatrix(value)
        if scaling_matrix.shape != self.lattice.matrix.shape:
            scaling_matrix = np.diagflat(scaling_matrix)
        self._scaling_matrix = scaling_matrix

        self.lattice = self.lattice.__class__(
            cell_matrix=self.scaling_matrix * self.lattice.matrix)

        tvecs = \
            np.asarray(
                np.asmatrix(supercell_lattice_points(self.scaling_matrix)) *
                self.lattice.matrix)

        basis = self.basis[:]
        self.basis = BasisAtoms()
        for atom in basis:
            for tvec in tvecs:
                xs, ys, zs = \
                    self.lattice.cartesian_to_fractional(atom.r + tvec)
                if self.wrap_coords:
                    xs, ys, zs = \
                        self.lattice.wrap_fractional_coordinate(
                            [xs, ys, zs])
                self.basis.append(BasisAtom(atom.element, lattice=self.lattice,
                                            xs=xs, ys=ys, zs=zs))

    def rotate(self, **kwargs):
        """Rotate crystal cell lattice, basis, and unit cell."""
        if self.lattice is not None:
            self.lattice.rotate(**kwargs)
        if self.basis is not None:
            self.basis.rotate(**kwargs)
        self.unit_cell.rotate(**kwargs)

    def translate(self, t, fix_anchor_points=True):
        """Translate crystal cell basis."""
        if not fix_anchor_points and self.lattice is not None:
            self.lattice.translate(t)
        if self.basis is not None:
            self.basis.translate(t, fix_anchor_points=fix_anchor_points)
        self.unit_cell.translate(t, fix_anchor_points=fix_anchor_points)

    def update_basis(self, element, index=None, step=None):
        """Update a crystal cell basis element."""
        if index is None:
            [self.unit_cell.basis.__setitem__(i, element)
             for i in range(len(self.unit_cell.basis))]
            [self.basis.__setitem__(i, element)
             for i in range(len(self.basis))]
        elif isinstance(index, int):
            if step is None:
                step = self.unit_cell.basis.Natoms
            [self.unit_cell.basis.__setitem__(i, element)
             for i in range(index, len(self.unit_cell.basis), step)]
            [self.basis.__setitem__(i, element)
             for i in range(index, len(self.basis), step)]
        elif isinstance(index, (list, np.ndarray)):
            [self.unit_cell.basis.__setitem__(i, element) for i in index]
            [self.basis.__setitem__(i, element) for i in index]

    def todict(self):
        try:
            return dict(lattice=self.lattice,
                        basis=self.basis.symbols.tolist(),
                        coords=self.basis.rs.tolist(),
                        wrap_coords=self.wrap_coords,
                        unit_cell=self.unit_cell,
                        scaling_matrix=self.scaling_matrix.tolist())
        except AttributeError:
            return dict(lattice=self.lattice, basis=None, coords=None,
                        wrap_coords=self.wrap_coords,
                        unit_cell=self.unit_cell,
                        scaling_matrix=self.scaling_matrix.tolist())


class SuperCell(CrystalCell):
    """Class representation of crystal structure supercell.

    Parameters
    ----------
    unit_cell : :class:`~sknano.core.crystallography.UnitCell`
    scaling_matrix : {:class:`~python:int`, :class:`~python:list`}

    """

    def __init__(self, unit_cell, scaling_matrix, wrap_coords=False):
        if not isinstance(unit_cell, UnitCell):
            raise ValueError('Expected a `UnitCell` for `unit_cell`.')
        if not isinstance(scaling_matrix,
                          (int, float, tuple, list, np.ndarray)):
            raise ValueError('Expected an `int` or `array_like` object of\n'
                             'integers for `scaling_matrix`')
        super().__init__(unit_cell=unit_cell, scaling_matrix=scaling_matrix,
                         wrap_coords=wrap_coords)
