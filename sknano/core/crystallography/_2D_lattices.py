# -*- coding: utf-8 -*-
"""
=============================================================================
2D crystal lattice classes (:mod:`sknano.core.crystallography._2D_lattices`)
=============================================================================

.. currentmodule:: sknano.core.crystallography._2D_lattices

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
# from builtins import super
from builtins import dict
from future import standard_library
standard_library.install_aliases()
# from future.utils import with_metaclass
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Point, Vector

from ._mixins import Direct2DLatticeMixin, Reciprocal2DLatticeMixin, \
    UnitCellMixin

__all__ = ['Crystal2DLattice', 'Reciprocal2DLattice']


class Crystal2DLattice(Reciprocal2DLatticeMixin, UnitCellMixin):
    """2D crystal lattice class."""

    def __init__(self, a=None, b=None, gamma=None,
                 a1=None, a2=None, cell_matrix=None,
                 orientation_matrix=None):

        self.offset = Point(np.zeros(2))

        if cell_matrix is not None:
            cell_matrix = np.asarray(cell_matrix)
            a1 = np.array(cell_matrix[0, :])
            a2 = np.array(cell_matrix[1, :])

        if a1 is not None and a2 is not None:
            a1 = Vector(a1)
            a2 = Vector(a2)
            a = a1.length
            b = a2.length
            gamma = np.degrees(a1.angle(a2))
            cell_matrix = np.matrix(np.vstack((np.asarray(a1),
                                               np.asarray(a2))))

        self.a = a
        self.b = b
        self.gamma = gamma

        if cell_matrix is not None:
            orientation_matrix = \
                cell_matrix.T * self.fractional_matrix

        if orientation_matrix is None:
            orientation_matrix = np.matrix(np.identity(2))

        self.orientation_matrix = orientation_matrix
        self.lattice_type = None

        self.fmtstr = "a={a!r}, b={b!r}, gamma={gamma!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def todict(self):
        """Return `dict` of `CrystalLattice` parameters."""
        return dict(a=self.a, b=self.b, gamma=self.gamma)

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
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        return np.around(np.cos(np.radians(self.gamma)), decimals=6)

    @property
    def sin_gamma(self):
        """:math:`\\sin\\gamma`"""
        return np.around(np.sin(np.radians(self.gamma)), decimals=6)

    @property
    def cell_area(self):
        """Unit cell area."""
        return self.a1.cross(self.a2)

    @property
    def ortho_matrix(self):
        """Transformation matrix to convert from fractional coordinates to \
            cartesian coordinates."""
        m11 = self.a
        m12 = self.b * self.cos_gamma
        m22 = self.b * self.sin_gamma

        return np.matrix([[m11, m12],
                          [0.0, m22]])

    @property
    def reciprocal_lattice(self):
        """Return `Crystal2DLattice` reciprocal lattice."""
        return Reciprocal2DLattice(cell_matrix=np.linalg.inv(self.cell_matrix))

    @classmethod
    def oblique(cls, a, b, gamma):
        """Generate an oblique 2D lattice with lattice parameters \
            :math:`a, b, \\gamma`."""
        return cls(a=a, b=b, gamma=gamma)

    @classmethod
    def rectangular(cls, a, b):
        """Generate a rectangular 2D lattice with lattice parameters \
            :math:`a, b`."""
        return cls(a=a, b=b, gamma=90)

    @classmethod
    def square(cls, a):
        """Generate a square 2D lattice with lattice parameter \
            :math:`a`."""
        return cls(a=a, b=a, gamma=90)

    @classmethod
    def hexagonal(cls, a):
        """Generate a hexagonal 2D lattice with lattice parameter \
            :math:`a`."""
        return cls(a=a, b=a, gamma=120)


class Reciprocal2DLattice(Direct2DLatticeMixin, UnitCellMixin):
    def __init__(self, a_star=None, b_star=None, gamma_star=None,
                 b1=None, b2=None, cell_matrix=None, orientation_matrix=None):

        if cell_matrix is not None:
            cell_matrix = np.asarray(cell_matrix)
            b1 = np.array(cell_matrix[0, :])
            b2 = np.array(cell_matrix[1, :])

        if b1 is not None and b2 is not None:
            b1 = Vector(b1)
            b2 = Vector(b2)
            a_star = b1.length
            b_star = b2.length
            gamma_star = np.degrees(b1.angle(b2))
            cell_matrix = np.matrix(np.vstack((np.asarray(b1),
                                               np.asarray(b2))))

        self.a_star = a_star
        self.b_star = b_star
        self.gamma_star = gamma_star

        if cell_matrix is not None:
            orientation_matrix = \
                cell_matrix.T * self.fractional_matrix

        if orientation_matrix is None:
            orientation_matrix = np.matrix(np.identity(2))

        self.orientation_matrix = orientation_matrix
        self.lattice_type = None

        self.fmtstr = "a_star={a_star!r}, b_star={b_star!r}, " + \
            "gamma_star={gamma_star!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def todict(self):
        """Return `dict` of `Reciprocal2DLattice` parameters."""
        return dict(a_star=self.a_star, b_star=self.b_star,
                    gamma_star=self.gamma_star)

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
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^*`"""
        return np.around(np.cos(np.radians(self.gamma_star)), decimals=6)

    @property
    def sin_gamma_star(self):
        """:math:`\\sin\\gamma^*`"""
        return np.around(np.sin(np.radians(self.gamma_star)), decimals=6)

    @property
    def cell_area(self):
        """Reciprocal unit cell area."""
        return self.b1.cross(self.b2)

    @property
    def ortho_matrix(self):
        """Transformation matrix to convert from fractional coordinates to \
            cartesian coordinates."""
        m11 = self.a_star
        m12 = self.b_star * self.cos_gamma_star
        m22 = self.b_star * self.sin_gamma_star
        return np.matrix([[m11, m12],
                          [0.0, m22]])

    @property
    def reciprocal_lattice(self):
        """`ReciprocalLattice` reciprocal lattice."""
        return Crystal2DLattice(cell_matrix=np.linalg.inv(self.cell_matrix))

    @classmethod
    def oblique(cls, a_star, b_star, gamma_star):
        """Generate an oblique 2D lattice with lattice parameters \
            :math:`a^*, b^*, \\gamma^*`."""
        return cls(a_star=a_star, b_star=b_star, gamma_star=gamma_star)

    @classmethod
    def rectangular(cls, a_star, b_star):
        """Generate a rectangular 2D lattice with lattice parameters \
            :math:`a^*, b^*`."""
        return cls(a_star=a_star, b_star=b_star, gamma_star=90)

    @classmethod
    def square(cls, a_star):
        """Generate a square 2D lattice with lattice parameter \
            :math:`a^*`."""
        return cls(a_star=a_star, b_star=a_star, gamma_star=90)

    @classmethod
    def hexagonal(cls, a_star):
        """Generate a hexagonal 2D lattice with lattice parameter \
            :math:`a^*`."""
        return cls(a_star=a_star, b_star=a_star, gamma_star=120)
