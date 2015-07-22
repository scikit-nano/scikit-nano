# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin crystallography classes (:mod:`sknano.core.crystallography._mixins`)
===============================================================================

.. currentmodule:: sknano.core.crystallography._mixins

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty
# from enum import Enum

from sknano.core.math import Point, Vector, zhat, rotation_matrix

import numpy as np

__all__ = ['DirectLatticeMixin', 'ReciprocalLatticeMixin', 'UnitCellMixin',
           'Direct2DLatticeMixin', 'Direct3DLatticeMixin',
           'Reciprocal2DLatticeMixin', 'Reciprocal3DLatticeMixin']


class Direct2DLatticeMixin:
    """Mixin class for computing the 2D direct lattice parameters from \
        2D reciprocal lattice parameters."""
    @property
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        return np.around(np.acos(np.radians(self.gamma)), decimals=10)

    @property
    def sin_gamma(self):
        """:math:`\\sin\\gamma`"""
        return np.sqrt(1 - self.cos_gamma ** 2)

    @property
    def a1(self):
        """2D lattice vector :math:`\\mathbf{a}_1=\\mathbf{a}`."""
        b2 = Vector()
        b2[:2] = self.b2
        return b2.cross(zhat) / self.cell_area

    @property
    def a2(self):
        """2D lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        b1 = Vector()
        b1[:2] = self.b1
        return zhat.cross(b1) / self.cell_area


class Reciprocal2DLatticeMixin:
    """Mixin class for computing the 2D reciprocal lattice parameters from \
        2D direct lattice parameters."""
    @property
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^*`"""
        return np.around(np.acos(np.radians(self.gamma_star)), decimals=10)

    @property
    def sin_gamma_star(self):
        """:math:`\\sin\\gamma^*`"""
        return np.sqrt(1 - self.cos_gamma_star ** 2)

    @property
    def b1(self):
        """2D reciprocal lattice vector :math:`\\mathbf{b}_1=\\mathbf{a}^{*}`.
        """
        a2 = Vector()
        a2[:2] = self.a2
        return a2.cross(zhat)[:2] / self.cell_area

    @property
    def b2(self):
        """2D reciprocal lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^{*}`.
        """
        a1 = Vector()
        a1[:2] = self.a1
        return zhat.cross(a1)[:2] / self.cell_area


class Direct3DLatticeMixin:
    """Mixin class for computing the 3D direct lattice parameters from \
        3D reciprocal lattice parameters."""
    @property
    def cos_alpha(self):
        """:math:`\\cos\\alpha`"""
        return np.around(
            (self.cos_beta_star * self.cos_gamma_star - self.cos_alpha_star) /
            (self.sin_beta_star * self.sin_gamma_star), decimals=10)

    @property
    def cos_beta(self):
        """:math:`\\cos\\beta`"""
        return np.around(
            (self.cos_gamma * self.cos_alpha - self.cos_beta) /
            (self.sin_gamma * self.sin_alpha), decimals=10)

    @property
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        return np.around(
            (self.cos_alpha_star * self.cos_beta_star - self.cos_gamma_star) /
            (self.sin_alpha_star * self.sin_beta_star), decimals=10)

    @property
    def sin_alpha(self):
        """:math:`\\sin\\alpha`"""
        return np.sqrt(1 - self.cos_alpha ** 2)

    @property
    def sin_beta(self):
        """:math:`\\sin\\beta`"""
        return np.sqrt(1 - self.cos_beta ** 2)

    @property
    def sin_gamma(self):
        """:math:`\\sin\\gamma`"""
        return np.sqrt(1 - self.cos_gamma ** 2)

    @property
    def a1(self):
        """3D lattice vector :math:`\\mathbf{a}_1=\\mathbf{a}`."""
        return self.b2.cross(self.b3) / self.cell_volume

    @property
    def a2(self):
        """3D lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        return self.b3.cross(self.b1) / self.cell_volume

    @property
    def a3(self):
        """3D lattice vector :math:`\\mathbf{a}_3=\\mathbf{c}`."""
        return self.b1.cross(self.b2) / self.cell_volume

DirectLatticeMixin = Direct3DLatticeMixin


class Reciprocal3DLatticeMixin:
    """Mixin class for computing the 3D reciprocal lattice parameters from \
        the 3D direct lattice parameters."""
    @property
    def cos_alpha_star(self):
        """:math:`\\cos\\alpha^*`"""
        return np.around((self.cos_beta * self.cos_gamma - self.cos_alpha) /
                         (self.sin_beta * self.sin_gamma), decimals=10)

    @property
    def cos_beta_star(self):
        """:math:`\\cos\\beta^*`"""
        return np.around((self.cos_gamma * self.cos_alpha - self.cos_beta) /
                         (self.sin_gamma * self.sin_alpha), decimals=10)

    @property
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^*`"""
        return np.around((self.cos_alpha * self.cos_beta - self.cos_gamma) /
                         (self.sin_alpha * self.sin_beta), decimals=10)

    @property
    def sin_alpha_star(self):
        """:math:`\\sin\\alpha^*`"""
        return np.sqrt(1 - self.cos_alpha_star ** 2)

    @property
    def sin_beta_star(self):
        """:math:`\\sin\\beta^*`"""
        return np.sqrt(1 - self.cos_beta_star ** 2)

    @property
    def sin_gamma_star(self):
        """:math:`\\sin\\gamma^*`"""
        return np.sqrt(1 - self.cos_gamma_star ** 2)

    @property
    def b1(self):
        """3D reciprocal lattice vector :math:`\\mathbf{b}_1=\\mathbf{a}^{*}`.
        """
        return self.a2.cross(self.a3) / self.cell_volume

    @property
    def b2(self):
        """3D reciprocal lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^{*}`.
        """
        return self.a3.cross(self.a1) / self.cell_volume

    @property
    def b3(self):
        """3D reciprocal lattice vector :math:`\\mathbf{b}_3=\\mathbf{c}^{*}`.
        """
        return self.a1.cross(self.a2) / self.cell_volume

ReciprocalLatticeMixin = Reciprocal3DLatticeMixin


class UnitCellMixin:
    """Mixin class for lattice unit cell."""
    @property
    def cell_matrix(self):
        """Matrix of lattice row vectors.

        Same as :attr:`Crystal2DLattice.ortho_matrix`\ .T or
        :attr:`Crystal3DLattice.ortho_matrix`\ .T.

        """
        return self.ortho_matrix.T * self.orientation_matrix.T

    @property
    def fractional_matrix(self):
        """Transformation matrix to convert from cartesian coordinates to \
            fractional coordinates."""
        return np.linalg.inv(self.ortho_matrix)

    @property
    def metric_tensor(self):
        """Metric tensor."""
        return self.ortho_matrix * self.ortho_matrix.T

    def fractional_to_cartesian(self, p):
        """Convert fractional coordinate to cartesian coordinate.

        Parameters
        ----------
        p : `Point`

        Returns
        -------
        `Point`

        """
        p = Point(p)
        c = self.orientation_matrix * self.ortho_matrix * \
            p.column_matrix + self.offset.column_matrix
        return p.__class__(c.A.flatten())

    def cartesian_to_fractional(self, p):
        """Convert cartesian coordinate to fractional coordinate.

        Parameters
        ----------
        p : `Point`

        Returns
        -------
        `Point`

        """
        p = Point(p)
        f = self.fractional_matrix * np.linalg.inv(self.orientation_matrix) * \
            (p - self.offset).column_matrix
        return p.__class__(f.A.flatten())

    def wrap_fractional_coordinate(self, p, epsilon=1e-8):
        """Wrap fractional coordinate to lie within unit cell.

        Parameters
        ----------
        p : `Point`

        Returns
        -------
        `Point`

        """
        p = Point(p)
        p = np.fmod(p, 1)
        p[np.where(p.__array__() < 0)] += 1
        p[np.where(p.__array__() > 1 - epsilon)] -= 1
        p[np.where(np.logical_or(
                   (p.__array__() > 1 - epsilon),
                   (p.__array__() < epsilon)))] = 0
        return p

    def wrap_cartesian_coordinate(self, p):
        """Wrap cartesian coordinate to lie within unit cell.

        Parameters
        ----------
        p : `Point`

        Returns
        -------
        `Point`

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
