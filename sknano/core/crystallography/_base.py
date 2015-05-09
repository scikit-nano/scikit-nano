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
from builtins import super
# from future.utils import with_metaclass

__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty
# from enum import Enum

from sknano.core.math import Point, Vector

import numpy as np

__all__ = ['CrystalLattice', 'CrystalStructure',
           'ReciprocalLatticeMixin']


class ReciprocalLatticeVectors:
    pass


class ReciprocalLatticeMixin:
    """Mixin class for computing the reciprocal lattice properties."""
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


class CrystalLattice(ReciprocalLatticeMixin):
    """Base class for crystal lattice systems."""

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

        self.pstr = "a={a!r}, b={b!r}, c={c!r}, " + \
            "alpha={alpha!r}, beta={beta!r}, gamma={gamma!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.pstr.format(**self.pdict()))

    def pdict(self):
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
    def cell_matrix(self):
        """Matrix of `CrystalLattice` lattice row vectors. \
            Same as :attr:`CrystalLattice.ortho_matrix`\ .T."""
        return self.ortho_matrix.T

    @property
    def cell_volume(self):
        """Unit cell volume."""
        return self.a * self.b * self.c * \
            np.sqrt(1 - self.cos_alpha ** 2 - self.cos_beta ** 2 -
                    self.cos_gamma ** 2 +
                    2 * self.cos_alpha * self.cos_beta * self.cos_gamma)

    @property
    def fractional_matrix(self):
        """Transformation matrix to convert from cartesian coordinates to \
            fractional coordinates."""
        return np.linalg.inv(self.ortho_matrix)

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

        p = Point(np.fmod(p, 1))
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

    # @property
    # def cell_vectors(self):
    #     pass

    # @property
    # def lattice_vectors(self):
    #     pass

    # @property
    # def lattice_type(self):
    #     return self._lattice_type

    # @lattice_type.setter
    # def lattice_type(self, value):
    #     self._lattice_type = value

    # @property
    # def origin(self):
    #     return self._origin

    # @origin.setter
    # def origin(self, value):
    #     self._origin = value

    # @property
    # def space_group(self):
    #     pass

    # def generate_cell_vectors(self):
    #     pass


class CrystalStructure(CrystalLattice):
    """Abstract base class for crystal structures."""

    def __init__(self, basis, **kwargs):
        super().__init__(**kwargs)
        self.basis = basis

        self.pstr = "basis={basis!r}"
        for k, v in kwargs.items():
            self.pstr += ", {}={{{key!s}!r}}".format(k, key=k)

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.pstr.format(**self.pdict()))

    @property
    def basis(self):
        """Crystal structure basis."""
        return self._basis

    @basis.setter
    def basis(self, value):
        self._basis = value

    @property
    def unit_cell(self):
        pass

    def pdict(self):
        """Return `dict` of `CrystalStructure` parameters."""
        super_pdict = super().pdict()
        super_pdict.update(dict(basis=self.basis))
        return super_pdict
