# -*- coding: utf-8 -*-
"""
===============================================================================
Base crystallography classes (:mod:`sknano.core.crystallography._base`)
===============================================================================

.. currentmodule:: sknano.core.crystallography._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
# from builtins import object
# from builtins import super
# from future.utils import with_metaclass

__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty
# from enum import Enum

from sknano.core.math import Point, Vector

import numpy as np

__all__ = ['DirectLatticeMixin', 'ReciprocalLatticeMixin', 'UnitCellMixin',
           'CrystalLattice', 'ReciprocalLattice', 'CrystalStructure']


class DirectLatticeMixin:
    """Mixin class for computing the direct lattice parameters from \
        reciprocal lattice parameters."""
    @property
    def cos_alpha(self):
        """:math:`\\cos\\alpha`"""
        return np.around(
            (self.cos_beta_star * self.cos_gamma_star - self.cos_alpha_star) /
            (self.sin_beta_star * self.sin_gamma_star), decimals=6)

    @property
    def cos_beta(self):
        """:math:`\\cos\\beta`"""
        return np.around(
            (self.cos_gamma * self.cos_alpha - self.cos_beta) /
            (self.sin_gamma * self.sin_alpha), decimals=6)

    @property
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        return np.around(
            (self.cos_alpha_star * self.cos_beta_star - self.cos_gamma_star) /
            (self.sin_alpha_star * self.sin_beta_star), decimals=6)

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
        """Lattice vector :math:`\\mathbf{a}_1=\\mathbf{a}`."""
        return self.b2.cross(self.b3) / self.cell_volume

    @property
    def a2(self):
        """Lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        return self.b3.cross(self.b1) / self.cell_volume

    @property
    def a3(self):
        """Lattice vector :math:`\\mathbf{a}_3=\\mathbf{c}`."""
        return self.b1.cross(self.b2) / self.cell_volume


class ReciprocalLatticeMixin:
    """Mixin class for computing the reciprocal lattice parameters from \
        the direct lattice parameters."""
    @property
    def cos_alpha_star(self):
        """:math:`\\cos\\alpha^*`"""
        return np.around((self.cos_beta * self.cos_gamma - self.cos_alpha) /
                         (self.sin_beta * self.sin_gamma), decimals=6)

    @property
    def cos_beta_star(self):
        """:math:`\\cos\\beta^*`"""
        return np.around((self.cos_gamma * self.cos_alpha - self.cos_beta) /
                         (self.sin_gamma * self.sin_alpha), decimals=6)

    @property
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^*`"""
        return np.around((self.cos_alpha * self.cos_beta - self.cos_gamma) /
                         (self.sin_alpha * self.sin_beta), decimals=6)

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
        """Reciprocal lattice vector :math:`\\mathbf{b}_1=\\mathbf{a}^{*}`."""
        return self.a2.cross(self.a3) / self.cell_volume

    @property
    def b2(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^{*}`."""
        return self.a3.cross(self.a1) / self.cell_volume

    @property
    def b3(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_3=\\mathbf{c}^{*}`."""
        return self.a1.cross(self.a2) / self.cell_volume


class UnitCellMixin:
    """Mixin class for lattice unit cell."""
    @property
    def cell_matrix(self):
        """Matrix of `CrystalLattice` lattice row vectors. \
            Same as :attr:`CrystalLattice.ortho_matrix`\ .T."""
        return self.ortho_matrix.T

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
        return Point(c.A.flatten())

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
        return Point(f.A.flatten())

    def wrap_fractional_coordinate(self, p, epsilon=1e-6):
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


class CrystalLattice(ReciprocalLatticeMixin, UnitCellMixin):
    """Crystal lattice class."""

    def __init__(self, a=None, b=None, c=None,
                 alpha=None, beta=None, gamma=None,
                 a1=None, a2=None, a3=None, cell_matrix=None,
                 orientation_matrix=None):

        self.offset = Point(np.zeros(3))

        if cell_matrix is not None:
            cell_matrix = np.asarray(cell_matrix)
            a1 = np.array(cell_matrix[0, :])
            a2 = np.array(cell_matrix[1, :])
            a3 = np.array(cell_matrix[2, :])

        if a1 is not None and a2 is not None and a3 is not None:
            a1 = Vector(a1)
            a2 = Vector(a2)
            a3 = Vector(a3)
            a = a1.length
            b = a2.length
            c = a3.length
            alpha = np.degrees(a2.angle(a3))
            beta = np.degrees(a3.angle(a1))
            gamma = np.degrees(a1.angle(a2))
            cell_matrix = np.matrix(np.vstack((np.asarray(a1),
                                               np.asarray(a2),
                                               np.asarray(a3))))

        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        if cell_matrix is not None:
            orientation_matrix = \
                cell_matrix.T * self.fractional_matrix

        if orientation_matrix is None:
            orientation_matrix = np.matrix(np.identity(3))

        self.orientation_matrix = orientation_matrix
        self.lattice_type = None

        self.fmtstr = "a={a!r}, b={b!r}, c={c!r}, " + \
            "alpha={alpha!r}, beta={beta!r}, gamma={gamma!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def todict(self):
        """Return `dict` of `CrystalLattice` parameters."""
        return dict(a=self.a, b=self.b, c=self.c,
                    alpha=self.alpha, beta=self.beta, gamma=self.gamma)

    @property
    def alpha(self):
        """Angle between lattice vectors :math:`\\mathbf{b}` and \
        :math:`\\mathbf{c}` in **degrees**."""
        return np.around(self._alpha, decimals=6)

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    @property
    def beta(self):
        """Angle between lattice vectors :math:`\\mathbf{c}` and \
        :math:`\\mathbf{a}` in **degrees**."""
        return np.around(self._beta, decimals=6)

    @beta.setter
    def beta(self, value):
        self._beta = value

    @property
    def gamma(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{b}` in **degrees**."""
        return np.around(self._gamma, decimals=6)

    @gamma.setter
    def gamma(self, value):
        self._gamma = value

    @property
    def a1(self):
        """Lattice vector :math:`\\mathbf{a}_1=\\mathbf{a}`."""
        return Vector(self.ortho_matrix[:, 0].A.flatten())

    @property
    def a2(self):
        """Lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        return Vector(self.ortho_matrix[:, 1].A.flatten())

    @property
    def a3(self):
        """Lattice vector :math:`\\mathbf{a}_3=\\mathbf{c}`."""
        return Vector(self.ortho_matrix[:, 2].A.flatten())

    @property
    def cos_alpha(self):
        """:math:`\\cos\\alpha`"""
        return np.around(np.cos(np.radians(self.alpha)), decimals=6)

    @property
    def cos_beta(self):
        """:math:`\\cos\\beta`"""
        return np.around(np.cos(np.radians(self.beta)), decimals=6)

    @property
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        return np.around(np.cos(np.radians(self.gamma)), decimals=6)

    @property
    def sin_alpha(self):
        """:math:`\\sin\\alpha`"""
        return np.around(np.sin(np.radians(self.alpha)), decimals=6)

    @property
    def sin_beta(self):
        """:math:`\\sin\\beta`"""
        return np.around(np.sin(np.radians(self.beta)), decimals=6)

    @property
    def sin_gamma(self):
        """:math:`\\sin\\gamma`"""
        return np.around(np.sin(np.radians(self.gamma)), decimals=6)

    @property
    def cell_volume(self):
        """Unit cell volume."""
        return self.a * self.b * self.c * \
            np.sqrt(1 - self.cos_alpha ** 2 - self.cos_beta ** 2 -
                    self.cos_gamma ** 2 +
                    2 * self.cos_alpha * self.cos_beta * self.cos_gamma)

    @property
    def ortho_matrix(self):
        """Transformation matrix to convert from fractional coordinates to \
            cartesian coordinates."""
        m11 = self.a
        m12 = self.b * self.cos_gamma
        m13 = self.c * self.cos_beta

        m22 = self.b * self.sin_gamma
        m23 = self.c * (self.cos_alpha - self.cos_beta * self.cos_gamma) / \
            self.sin_gamma

        m33 = self.c * self.sin_alpha * self.sin_beta * self.sin_gamma_star / \
            self.sin_gamma

        return np.matrix([[m11, m12, m13],
                          [0.0, m22, m23],
                          [0.0, 0.0, m33]])

    @property
    def reciprocal_lattice(self):
        """Return `CrystalLattice` reciprocal lattice."""
        return ReciprocalLattice(cell_matrix=np.linalg.inv(self.cell_matrix))


class ReciprocalLattice(DirectLatticeMixin, UnitCellMixin):
    def __init__(self, a_star=None, b_star=None, c_star=None,
                 alpha_star=None, beta_star=None, gamma_star=None,
                 b1=None, b2=None, b3=None,
                 cell_matrix=None, orientation_matrix=None):

        if cell_matrix is not None:
            cell_matrix = np.asarray(cell_matrix)
            b1 = np.array(cell_matrix[0, :])
            b2 = np.array(cell_matrix[1, :])
            b3 = np.array(cell_matrix[2, :])

        if b1 is not None and b2 is not None and b3 is not None:
            b1 = Vector(b1)
            b2 = Vector(b2)
            b3 = Vector(b3)
            a_star = b1.length
            b_star = b2.length
            c_star = b3.length
            alpha_star = np.degrees(b2.angle(b3))
            beta_star = np.degrees(b3.angle(b1))
            gamma_star = np.degrees(b1.angle(b2))
            cell_matrix = np.matrix(np.vstack((np.asarray(b1),
                                               np.asarray(b2),
                                               np.asarray(b3))))

        self.a_star = a_star
        self.b_star = b_star
        self.c_star = c_star
        self.alpha_star = alpha_star
        self.beta_star = beta_star
        self.gamma_star = gamma_star

        if cell_matrix is not None:
            orientation_matrix = \
                cell_matrix.T * self.fractional_matrix

        if orientation_matrix is None:
            orientation_matrix = np.matrix(np.identity(3))

        self.orientation_matrix = orientation_matrix
        self.lattice_type = None

        self.fmtstr = "a_star={a_star!r}, b_star={b_star!r}, " + \
            "c_star={c_star!r}, alpha_star={alpha_star!r}, " + \
            "beta_star={beta_star!r}, gamma_star={gamma_star!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def todict(self):
        """Return `dict` of `ReciprocalLattice` parameters."""
        return dict(a_star=self.a_star, b_star=self.b_star,
                    c_star=self.c_star, alpha_star=self.alpha_star,
                    beta_star=self.beta_star, gamma_star=self.gamma_star)

    @property
    def alpha_star(self):
        """Angle between lattice vectors :math:`\\mathbf{b}^{\\ast}` and \
        :math:`\\mathbf{c}^{\\ast}` in **degrees**."""
        return np.around(self._alpha_star, decimals=6)

    @alpha_star.setter
    def alpha_star(self, value):
        self._alpha_star = value

    @property
    def beta_star(self):
        """Angle between lattice vectors :math:`\\mathbf{c}` and \
        :math:`\\mathbf{a}` in **degrees**."""
        return np.around(self._beta_star, decimals=6)

    @beta_star.setter
    def beta_star(self, value):
        self._beta_star = value

    @property
    def gamma_star(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{b}` in **degrees**."""
        return np.around(self._gamma_star, decimals=6)

    @gamma_star.setter
    def gamma_star(self, value):
        self._gamma_star = value

    @property
    def b1(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_1=\\mathbf{a}^*`."""
        return Vector(self.ortho_matrix[:, 0].A.flatten())

    @property
    def b2(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^*`."""
        return Vector(self.ortho_matrix[:, 1].A.flatten())

    @property
    def b3(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_3=\\mathbf{c}^*`."""
        return Vector(self.ortho_matrix[:, 2].A.flatten())

    @property
    def cos_alpha_star(self):
        """:math:`\\cos\\alpha^*`"""
        return np.around(np.cos(np.radians(self.alpha_star)), decimals=6)

    @property
    def cos_beta_star(self):
        """:math:`\\cos\\beta^*`"""
        return np.around(np.cos(np.radians(self.beta_star)), decimals=6)

    @property
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^*`"""
        return np.around(np.cos(np.radians(self.gamma_star)), decimals=6)

    @property
    def sin_alpha_star(self):
        """:math:`\\sin\\alpha^*`"""
        return np.around(np.sin(np.radians(self.alpha_star)), decimals=6)

    @property
    def sin_beta_star(self):
        """:math:`\\sin\\beta^*`"""
        return np.around(np.sin(np.radians(self.beta_star)), decimals=6)

    @property
    def sin_gamma_star(self):
        """:math:`\\sin\\gamma^*`"""
        return np.around(np.sin(np.radians(self.gamma_star)), decimals=6)

    @property
    def cell_volume(self):
        """Reciprocal unit cell volume."""
        return self.a_star * self.b_star * self.c_star * \
            np.sqrt(1 - self.cos_alpha_star ** 2 - self.cos_beta_star ** 2 -
                    self.cos_gamma_star ** 2 + 2 * self.cos_alpha_star *
                    self.cos_beta_star * self.cos_gamma_star)

    @property
    def ortho_matrix(self):
        """Transformation matrix to convert from fractional coordinates to \
            cartesian coordinates."""
        m11 = self.a_star
        m12 = self.b_star * self.cos_gamma_star
        m13 = self.c_star * self.cos_beta_star

        m22 = self.b_star * self.sin_gamma_star
        m23 = self.c_star * (self.cos_alpha_star - self.cos_beta_star *
                             self.cos_gamma_star) / self.sin_gamma_star

        m33 = self.c_star * self.sin_alpha_star * self.sin_beta_star * \
            self.sin_gamma / self.sin_gamma_star

        return np.matrix([[m11, m12, m13],
                          [0.0, m22, m23],
                          [0.0, 0.0, m33]])

    @property
    def reciprocal_lattice(self):
        """`ReciprocalLattice` reciprocal lattice."""
        return CrystalLattice(cell_matrix=np.linalg.inv(self.cell_matrix))


class CrystalStructure:
    """Abstract base class for crystal structures."""

    def __init__(self, lattice, basis, coords=None, cartesian=False):

        self.lattice = lattice
        self.basis = basis

        self.fmtstr = "{lattice!r}, {basis!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    @property
    def basis(self):
        """Crystal structure basis."""
        return self._basis

    @basis.setter
    def basis(self, value):
        self._basis = value

    @property
    def lattice(self):
        return self._lattice

    @lattice.setter
    def lattice(self, value):
        if not isinstance(value, CrystalLattice):
            value = CrystalLattice(cell_matrix=value)
        self._lattice = value

    @property
    def unit_cell(self):
        pass

    def todict(self):
        """Return `dict` of `CrystalStructure` parameters."""
        return dict(lattice=self.lattice, basis=self.basis)
