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
from sknano.core.atoms import BasisAtom as Atom, BasisAtoms as Atoms
from sknano.core.math import rotation_matrix
from ._extras import supercell_lattice_points

__all__ = ['CrystalCellMixin', 'CrystalCell', 'UnitCell', 'SuperCell']


class CrystalCellMixin:
    """Mixin class for crystal lattice cell."""

    def fractional_to_cartesian(self, fcoords):
        """Convert fractional coordinate to cartesian coordinate.

        Parameters
        ----------
        fcoords : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        ccoords = self.orientation_matrix * self.ortho_matrix * \
            np.asmatrix(fcoords).T + self.offset.column_matrix
        try:
            return ccoords.T.A.reshape((3, ))
        except ValueError:
            return ccoords.T.A.reshape((len(fcoords), 3))

    def cartesian_to_fractional(self, ccoords):
        """Convert cartesian coordinate to fractional coordinate.

        Parameters
        ----------
        ccoords : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        fcoords = np.linalg.inv(self.ortho_matrix) * \
            np.linalg.inv(self.orientation_matrix) * \
            (np.asmatrix(ccoords).T - self.offset.column_matrix)
        try:
            return fcoords.T.A.reshape((3, ))
        except ValueError:
            return fcoords.T.A.reshape((len(ccoords), 3))

    def wrap_fractional_coordinate(self, p, epsilon=1e-8):
        """Wrap fractional coordinate to lie within unit cell.

        Parameters
        ----------
        p : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        p = np.asarray(p)
        p = np.fmod(p, 1)
        p[np.where(p < 0)] += 1
        p[np.where(p > 1 - epsilon)] -= 1
        p[np.where(np.logical_or((p > 1 - epsilon), (p < epsilon)))] = 0
        return p

    def wrap_cartesian_coordinate(self, p):
        """Wrap cartesian coordinate to lie within unit cell.

        Parameters
        ----------
        p : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        return self.fractional_to_cartesian(
            self.wrap_fractional_coordinate(
                self.cartesian_to_fractional(p)))

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None, degrees=False,
               transform_matrix=None, verbose=False, **kwargs):
        """Rotate unit cell.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        See Also
        --------
        core.math.rotate

        """
        if self.nd == 2:
            axis = 'z'
        if transform_matrix is None:
            transform_matrix = \
                rotation_matrix(angle=angle, axis=axis,
                                anchor_point=anchor_point,
                                rot_point=rot_point,
                                from_vector=from_vector,
                                to_vector=to_vector, degrees=degrees,
                                verbose=verbose, **kwargs)

            # transform_matrix = \
            #     transformation_matrix(angle=angle, axis=axis,
            #                           anchor_point=anchor_point,
            #                           rot_point=rot_point,
            #                           from_vector=from_vector,
            #                           to_vector=to_vector, degrees=degrees,
            #                           verbose=verbose, **kwargs)

        self.orientation_matrix = \
            np.dot(transform_matrix, self.orientation_matrix)

    def translate(self, t):
        """Translate unit cell.

        Parameters
        ----------
        t : :class:`Vector`

        See Also
        --------
        core.math.translate

        """
        self.offset.translate(t)


@total_ordering
class UnitCell(BaseClass):
    """Base class for abstract representations of crystallographic unit cells.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.LatticeBase` sub-class
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional

    """

    def __init__(self, lattice=None, basis=None, coords=None, cartesian=False):

        super().__init__()

        if basis is None:
            basis = Atoms()
        elif basis is not None:
            basis = Atoms(basis)
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
        self.fmtstr = "{lattice!r}, {basis!r}, {coords!r}, cartesian=False"

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
    #         value = Atoms()
    #     elif value is not None:
    #         value = Atoms(value, lattice=lattice)
    #         if coords is not None:
    #             for atom, pos in zip(basis, coords):
    #                 atom.lattice = lattice
    #                 if not cartesian:
    #                     atom.rs = pos
    #                 else:
    #                     atom.rs = lattice.cartesian_to_fractional(pos)

    def rotate(self, **kwargs):
        """Rotate lattice vectors and basis."""
        self.lattice.rotate(**kwargs)
        self.basis.rotate(**kwargs)

    def todict(self):
        """Return `dict` of `UnitCell` parameters."""
        return dict(lattice=self.lattice, basis=self.basis.symbols.tolist(),
                    coords=self.basis.rs.tolist())


@total_ordering
class CrystalCell(BaseClass):
    """Class representation of crystal structure cell.

    Parameters
    ----------
    unit_cell : :class:`~sknano.core.crystallography.UnitCell`
    scaling_matrix : {:class:`~python:int`, :class:`~python:list`}

    """

    def __init__(self, unit_cell=None, scaling_matrix=None, wrap_coords=False):
        super().__init__()

        # These attributes may be reset in the `@scaling_matrix.setter` method
        # and so they need to be initialized *before* setting
        # `self.scaling_matrix`.
        self._lattice = None
        self._basis = None

        if unit_cell is not None:
            self.unit_cell = unit_cell

        self.scaling_matrix = scaling_matrix
        self.wrap_coords = wrap_coords
        self.fmtstr = \
            "unit_cell={unit_cell!r}, scaling_matrix={scaling_matrix!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['unit_cell', 'scaling_matrix'])
        return attrs

    def __eq__(self, other):
        if all([mat is not None for mat in
                (self.scaling_matrix, other.scaling_matrix)]) and \
                self.scaling_matrix.shape == other.scaling_matrix.shape:
            return self is other or \
                (self.unit_cell == other.unit_cell and
                 np.allclose(self.scaling_matrix, other.scaling_matrix))
        else:
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
        if name != 'basis' and self.basis is not None:
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
        return self._basis if self._basis is not None else \
            self.unit_cell.basis

    @basis.setter
    def basis(self, value):
        self._basis = value
        self.unit_cell.basis[:] = self.basis[:self.unit_cell.basis.Natoms]

    @property
    def lattice(self):
        return self._lattice if self._lattice is not None else \
            self.unit_cell.lattice

    @property
    def unit_cell(self):
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, value):
        if not isinstance(value, UnitCell):
            raise ValueError('Expected a `UnitCell` object')
        self._unit_cell = value

    @property
    def scaling_matrix(self):
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
                value = self.unit_cell.lattice.nd * [int(value)]

        scaling_matrix = np.asmatrix(value, dtype=int)
        if scaling_matrix.shape != self.unit_cell.lattice.matrix.shape:
            scaling_matrix = np.diagflat(scaling_matrix)
        self._scaling_matrix = scaling_matrix
        self._lattice = self.unit_cell.lattice.__class__(
            cell_matrix=self.scaling_matrix * self.unit_cell.lattice.matrix)

        tvecs = np.asmatrix(supercell_lattice_points(self.scaling_matrix)) * \
            self.lattice.matrix

        self._basis = Atoms()
        for atom in self.unit_cell.basis:
            for tvec in tvecs:
                xs, ys, zs = \
                    self.lattice.cartesian_to_fractional(atom.r + tvec)
                if self.wrap_coords:
                    xs, ys, zs = \
                        self.lattice.wrap_fractional_coordinate([xs, ys, zs])
                self.basis.append(Atom(atom.element, lattice=self.lattice,
                                       xs=xs, ys=ys, zs=zs))

    def rotate(self, **kwargs):
        """Rotate lattice vectors and basis."""
        self.lattice.rotate(**kwargs)
        self.basis.rotate(**kwargs)
        self.unit_cell.rotate(**kwargs)

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
        return dict(unit_cell=self.unit_cell,
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

        self.fmtstr = ', '.join((self.unit_cell.fmtstr, super().fmtstr))

    def todict(self):
        """Return `dict` of `SuperCell` parameters."""
        return dict(lattice=self.lattice, basis=self.basis.symbols.tolist(),
                    coords=self.basis.rs.tolist(), unit_cell=self.unit_cell,
                    scaling_matrix=self.scaling_matrix.tolist())
