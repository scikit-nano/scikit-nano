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

import numpy as np

from sknano.core import BaseClass
from sknano.core.atoms import BasisAtoms as Atoms
from sknano.core.math import rotation_matrix

__all__ = ['UnitCellMixin', 'UnitCell']


class UnitCellMixin:
    """Mixin class for lattice unit cell."""
    @property
    def cell_matrix(self):
        """Matrix of lattice row vectors.

        Same as :attr:`Crystal2DLattice.ortho_matrix`\ .T or
        :attr:`Crystal3DLattice.ortho_matrix`\ .T.

        """
        return (self.orientation_matrix * self.ortho_matrix).T

    @property
    def matrix(self):
        """Alias for \
            :attr:`~sknano.core.crystallography.UnitCellMixin.cell_matrix`."""
        return self.cell_matrix

    @property
    def fractional_matrix(self):
        """Transformation matrix to convert from cartesian coordinates to \
            fractional coordinates."""
        return np.linalg.inv(self.ortho_matrix)

    @property
    def metric_tensor(self):
        """Metric tensor."""
        return self.ortho_matrix.T * self.ortho_matrix

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
