# -*- coding: utf-8 -*-
"""
=============================================================================
2D crystal lattice classes (:mod:`sknano.core.crystallography._2D_lattices`)
=============================================================================

.. currentmodule:: sknano.core.crystallography._2D_lattices

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Vector

from ._base import LatticeBase, ReciprocalLatticeBase
from ._mixins import Direct2DLatticeMixin, Reciprocal2DLatticeMixin, \
    UnitCellMixin

__all__ = ['Crystal2DLattice', 'Reciprocal2DLattice']


class Crystal2DLattice(LatticeBase, Reciprocal2DLatticeMixin, UnitCellMixin):
    """2D crystal lattice class.

    Parameters
    ----------
    a, b : float
    gamma : float
    a1, a2 : array_like
    cell_matrix : array_like
    orientation_matrix : array_like, optional

    """

    def __init__(self, a=None, b=None, gamma=None,
                 a1=None, a2=None, cell_matrix=None,
                 orientation_matrix=None):

        if cell_matrix is not None:
            cell_matrix = np.asarray(cell_matrix)
            a1 = np.array(cell_matrix[0, :])
            a2 = np.array(cell_matrix[1, :])

        if a1 is not None and a2 is not None:
            a1 = Vector(a1, nd=3)
            a2 = Vector(a2, nd=3)
            a = a1.length
            b = a2.length
            gamma = np.degrees(a1.angle(a2))
            cell_matrix = np.matrix(np.vstack((np.asarray(a1),
                                               np.asarray(a2),
                                               np.asarray([0, 0, 1]))))

        self._a = a
        self._b = b
        self._gamma = gamma

        if None not in (a, b, gamma):
            self._update_ortho_matrix()

        super().__init__(nd=2, cell_matrix=cell_matrix,
                         orientation_matrix=orientation_matrix)

        self.fmtstr = "a={a!r}, b={b!r}, gamma={gamma!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['a', 'b', 'gamma'])
        return attrs

    def todict(self):
        """Return `dict` of `Crystal2DLattice` parameters."""
        return dict(a=self.a, b=self.b, gamma=self.gamma)

    @property
    def a(self):
        """Length of lattice vector :math:`\\mathbf{a}`."""
        return np.around(self._a, decimals=10)

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def b(self):
        """Length of lattice_vector :math:`\\mathbf{b}`."""
        return np.around(self._b, decimals=10)

    @b.setter
    def b(self, value):
        self._b = value

    @property
    def gamma(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{b}` in **degrees**."""
        try:
            return np.around(float(self._gamma), decimals=10)
        except TypeError:
            return None

    @property
    def a1(self):
        """Lattice vector :math:`\\mathbf{a}_1=\\mathbf{a}`."""
        return Vector(self.cell_matrix[0, :].A.flatten())[:2]

    # @a1.setter
    # def a1(self, value):
    #     if not isinstance(value, Vector):
    #         value = Vector(value)
    #     self.cell_matrix[:, 0] = np.asarray(value)

    @property
    def a2(self):
        """Lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        return Vector(self.cell_matrix[1, :].A.flatten())[:2]

    @property
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        return np.around(np.cos(np.radians(self.gamma)), decimals=10)

    @property
    def sin_gamma(self):
        """:math:`\\sin\\gamma`"""
        return np.around(np.sin(np.radians(self.gamma)), decimals=10)

    @property
    def cell_area(self):
        """Alias for :attr:`~Crystal2DLattice.area`."""
        return self.area

    @property
    def area(self):
        """Unit cell area."""
        return np.abs(self.a1.cross(self.a2))

    @property
    def ortho_matrix(self):
        """Transformation matrix to convert from fractional coordinates to \
            cartesian coordinates."""
        return self._ortho_matrix

    def _update_ortho_matrix(self):
        m11 = self.a
        m12 = self.b * self.cos_gamma
        m22 = self.b * self.sin_gamma

        # self._ortho_matrix = np.matrix([[m11, m12],
        #                                 [0.0, m22]])
        self._ortho_matrix = np.matrix([[m11, m12, 0.0],
                                        [0.0, m22, 0.0],
                                        [0.0, 0.0, 1.0]])

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


class Reciprocal2DLattice(ReciprocalLatticeBase, Direct2DLatticeMixin,
                          UnitCellMixin):
    """2D reciprocal lattice class.

    Parameters
    ----------
    a_star, b_star : float
    gamma_star : float
    b1, b2 : array_like
    cell_matrix : array_like
    orientation_matrix : array_like

    """

    def __init__(self, a_star=None, b_star=None, gamma_star=None,
                 b1=None, b2=None, cell_matrix=None, orientation_matrix=None):

        direct_lattice = \
            Crystal2DLattice(a=a_star, b=b_star, gamma=gamma_star,
                             a1=b1, a2=b2, cell_matrix=cell_matrix,
                             orientation_matrix=orientation_matrix)
        super().__init__(direct_lattice, nd=2)

        self.fmtstr = "a_star={a_star!r}, b_star={b_star!r}, " + \
            "gamma_star={gamma_star!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['a_star', 'b_star', 'gamma_star'])
        return attrs

    def todict(self):
        """Return `dict` of `Reciprocal2DLattice` parameters."""
        return dict(a_star=self.a_star, b_star=self.b_star,
                    gamma_star=self.gamma_star)

    @property
    def a_star(self):
        """Length of reciprocal lattice vector :math:`\\mathbf{a^*}`."""
        return self._direct_lattice.a

    @a_star.setter
    def a_star(self, value):
        self._direct_lattice.a = value

    @property
    def b_star(self):
        """Length of reciprocal lattice_vector :math:`\\mathbf{b^*}`."""
        return self._direct_lattice.b

    @b_star.setter
    def b_star(self, value):
        self._direct_lattice.b = value

    @property
    def gamma_star(self):
        """Angle between reciprocal lattice vectors \
            :math:`\\mathbf{a}^{\\ast}` and :math:`\\mathbf{b}^{\\ast}` in \
            **degrees**."""
        try:
            return np.around(self._direct_lattice.gamma, decimals=10)
        except TypeError:
            return None

    @property
    def b1(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_1=\\mathbf{a}^*`."""
        return self._direct_lattice.a1

    @property
    def b2(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^*`."""
        return self._direct_lattice.a2

    @property
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^*`"""
        return self._direct_lattice.cos_gamma

    @property
    def sin_gamma_star(self):
        """:math:`\\sin\\gamma^*`"""
        return self._direct_lattice.sin_gamma

    @property
    def reciprocal_lattice(self):
        """Reciprocal lattice of this `Reciprocal2DLattice`."""
        return Crystal2DLattice(cell_matrix=np.linalg.inv(self.cell_matrix))

    @classmethod
    def oblique(cls, a_star, b_star, gamma_star):
        """Generate an oblique 2D reciprocal lattice with lattice parameters \
            :math:`a^*, b^*, \\gamma^*`."""
        return cls(a_star=a_star, b_star=b_star, gamma_star=gamma_star)

    @classmethod
    def rectangular(cls, a_star, b_star):
        """Generate a rectangular 2D reciprocal lattice with lattice \
            parameters :math:`a^*, b^*`."""
        return cls(a_star=a_star, b_star=b_star, gamma_star=90)

    @classmethod
    def square(cls, a_star):
        """Generate a square 2D reciprocal lattice with lattice parameter \
            :math:`a^*`."""
        return cls(a_star=a_star, b_star=a_star, gamma_star=90)

    @classmethod
    def hexagonal(cls, a_star):
        """Generate a hexagonal 2D reciprocal lattice with lattice parameter \
            :math:`a^*`."""
        return cls(a_star=a_star, b_star=a_star, gamma_star=120)
