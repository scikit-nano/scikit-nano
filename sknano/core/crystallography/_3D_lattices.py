# -*- coding: utf-8 -*-
"""
=============================================================================
3D crystal lattice classes (:mod:`sknano.core.crystallography._3D_lattices`)
=============================================================================

.. currentmodule:: sknano.core.crystallography._3D_lattices

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Vector

from ._base import LatticeBase, ReciprocalLatticeBase
from ._mixins import DirectLatticeMixin, ReciprocalLatticeMixin, UnitCellMixin

__all__ = ['CrystalLattice', 'ReciprocalLattice',
           'Crystal3DLattice', 'Reciprocal3DLattice']
# __all__ += ['BravaisLattice', 'Bravais3DLattice',
#            'SimpleCubicLattice', 'BodyCenteredCubicLattice',
#            'FaceCenteredCubicLattice']


class Crystal3DLattice(LatticeBase, ReciprocalLatticeMixin, UnitCellMixin):
    """3D crystal lattice class.

    Parameters
    ----------
    a, b, c : float
    alpha, beta, gamma : float
    a1, a2, a3 : array_like
    cell_matrix : array_like
    orientation_matrix : array_like, optional

    """

    def __init__(self, a=None, b=None, c=None,
                 alpha=None, beta=None, gamma=None,
                 a1=None, a2=None, a3=None, cell_matrix=None,
                 orientation_matrix=None):

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

        self._a = a
        self._b = b
        self._c = c
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma

        if None not in (a, b, c, alpha, beta, gamma):
            self._update_ortho_matrix()

        super().__init__(nd=3, cell_matrix=cell_matrix,
                         orientation_matrix=orientation_matrix)

        self.fmtstr = "a={a!r}, b={b!r}, c={c!r}, " + \
            "alpha={alpha!r}, beta={beta!r}, gamma={gamma!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['a', 'b', 'c', 'alpha', 'beta', 'gamma'])
        return attrs

    def todict(self):
        """Return `dict` of `Crystal3DLattice` parameters."""
        return dict(a=self.a, b=self.b, c=self.c,
                    alpha=self.alpha, beta=self.beta, gamma=self.gamma)

    @property
    def a(self):
        """Length of lattice vector :math:`\\mathbf{a}`."""
        return np.around(self._a, decimals=10)

    @a.setter
    def a(self, value):
        self._a = value
        self._update_ortho_matrix()

    @property
    def b(self):
        """Length of lattice_vector :math:`\\mathbf{b}`."""
        return np.around(self._b, decimals=10)

    @b.setter
    def b(self, value):
        self._b = value
        self._update_ortho_matrix()

    @property
    def c(self):
        """Length of lattice vector :math:`\\mathbf{c}`."""
        return np.around(self._c, decimals=10)

    @c.setter
    def c(self, value):
        self._c = value
        self._update_ortho_matrix()

    @property
    def alpha(self):
        """Angle between lattice vectors :math:`\\mathbf{b}` and \
        :math:`\\mathbf{c}` in **degrees**."""
        return np.around(self._alpha, decimals=10)

    @property
    def beta(self):
        """Angle between lattice vectors :math:`\\mathbf{c}` and \
        :math:`\\mathbf{a}` in **degrees**."""
        return np.around(self._beta, decimals=10)

    @property
    def gamma(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{b}` in **degrees**."""
        return np.around(self._gamma, decimals=10)

    @property
    def a1(self):
        """Lattice vector :math:`\\mathbf{a}_1=\\mathbf{a}`."""
        return Vector(self.cell_matrix[0, :].A.flatten())

    @property
    def a2(self):
        """Lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        return Vector(self.cell_matrix[1, :].A.flatten())

    @property
    def a3(self):
        """Lattice vector :math:`\\mathbf{a}_3=\\mathbf{c}`."""
        return Vector(self.cell_matrix[2, :].A.flatten())

    @property
    def cos_alpha(self):
        """:math:`\\cos\\alpha`"""
        return np.around(np.cos(np.radians(self.alpha)), decimals=10)

    @property
    def cos_beta(self):
        """:math:`\\cos\\beta`"""
        return np.around(np.cos(np.radians(self.beta)), decimals=10)

    @property
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        return np.around(np.cos(np.radians(self.gamma)), decimals=10)

    @property
    def sin_alpha(self):
        """:math:`\\sin\\alpha`"""
        return np.around(np.sin(np.radians(self.alpha)), decimals=10)

    @property
    def sin_beta(self):
        """:math:`\\sin\\beta`"""
        return np.around(np.sin(np.radians(self.beta)), decimals=10)

    @property
    def sin_gamma(self):
        """:math:`\\sin\\gamma`"""
        return np.around(np.sin(np.radians(self.gamma)), decimals=10)

    @property
    def cell_volume(self):
        """Alias for :attr:`~Crystal3DLattice.volume`."""
        return self.volume

    @property
    def volume(self):
        """Unit cell volume."""
        return self.a * self.b * self.c * \
            np.sqrt(1 - self.cos_alpha ** 2 - self.cos_beta ** 2 -
                    self.cos_gamma ** 2 +
                    2 * self.cos_alpha * self.cos_beta * self.cos_gamma)

    @property
    def ortho_matrix(self):
        """Transformation matrix to convert from fractional coordinates to \
            cartesian coordinates."""
        return self._ortho_matrix

    def _update_ortho_matrix(self):
        m11 = self.a
        m12 = self.b * self.cos_gamma
        m13 = self.c * self.cos_beta

        m22 = self.b * self.sin_gamma
        m23 = self.c * (self.cos_alpha - self.cos_beta * self.cos_gamma) / \
            self.sin_gamma

        m33 = self.c * self.sin_alpha * self.sin_beta * self.sin_gamma_star / \
            self.sin_gamma

        self._ortho_matrix = np.matrix([[m11, m12, m13],
                                        [0.0, m22, m23],
                                        [0.0, 0.0, m33]])

    @property
    def reciprocal_lattice(self):
        """Reciprocal lattice of this `Crystal3DLattice`."""
        return Reciprocal3DLattice(cell_matrix=np.linalg.inv(self.cell_matrix))

    @classmethod
    def triclinic(cls, a, b, c, alpha, beta, gamma):
        """Generate a triclinic 3D lattice with lattice parameters \
            :math:`a, b, c, \\alpha, \\beta, \\gamma`."""
        return cls(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)

    @classmethod
    def monoclinic(cls, a, b, c, beta):
        """Generate a monoclinic 3D lattice with lattice parameters \
            :math:`a, b, c, \\beta`."""
        return cls(a=a, b=b, c=c, alpha=90, beta=beta, gamma=90)

    @classmethod
    def orthorhombic(cls, a, b, c):
        """Generate an orthorhombic 3D lattice with lattice parameters \
            :math:`a, b, c`."""
        return cls(a=a, b=b, c=c, alpha=90, beta=90, gamma=90)

    @classmethod
    def tetragonal(cls, a, c):
        """Generate a tetragonal 3D lattice with lattice parameters \
            :math:`a, c`."""
        return cls(a=a, b=a, c=c, alpha=90, beta=90, gamma=90)

    @classmethod
    def hexagonal(cls, a, c):
        """Generate a hexagonal 3D lattice with lattice parameters \
            :math:`a, c`."""
        return cls(a=a, b=a, c=c, alpha=90, beta=90, gamma=120)

    @classmethod
    def cubic(cls, a):
        """Generate a cubic 3D lattice with lattice parameter \
            :math:`a`."""
        return cls(a=a, b=a, c=a, alpha=90, beta=90, gamma=90)

CrystalLattice = Crystal3DLattice


class Reciprocal3DLattice(ReciprocalLatticeBase, DirectLatticeMixin,
                          UnitCellMixin):
    """3D reciprocal lattice class.

    Parameters
    ----------
    a_star, b_star, c_star : float
    alpha_star, beta_star, gamma_star : float
    b1, b2, b3 : array_like
    cell_matrix : array_like
    orientation_matrix : array_like

    """
    def __init__(self, a_star=None, b_star=None, c_star=None,
                 alpha_star=None, beta_star=None, gamma_star=None,
                 b1=None, b2=None, b3=None,
                 cell_matrix=None, orientation_matrix=None):

        direct_lattice = \
            Crystal3DLattice(a=a_star, b=b_star, c=c_star, alpha=alpha_star,
                             beta=beta_star, gamma=gamma_star,
                             a1=b1, a2=b2, a3=b3, cell_matrix=cell_matrix,
                             orientation_matrix=orientation_matrix)

        super().__init__(direct_lattice, nd=3)

        self.fmtstr = "a_star={a_star!r}, b_star={b_star!r}, " + \
            "c_star={c_star!r}, alpha_star={alpha_star!r}, " + \
            "beta_star={beta_star!r}, gamma_star={gamma_star!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['a_star', 'b_star', 'c_star',
                      'alpha_star', 'beta_star', 'gamma_star'])
        return attrs

    def todict(self):
        """Return `dict` of `Reciprocal3DLattice` parameters."""
        return dict(a_star=self.a_star, b_star=self.b_star,
                    c_star=self.c_star, alpha_star=self.alpha_star,
                    beta_star=self.beta_star, gamma_star=self.gamma_star)

    @property
    def a_star(self):
        """Length of lattice vector :math:`\\mathbf{a^*}`."""
        return self._direct_lattice.a

    @a_star.setter
    def a_star(self, value):
        self._direct_lattice.a = value

    @property
    def b_star(self):
        """Length of lattice_vector :math:`\\mathbf{b^*}`."""
        return self._direct_lattice.b

    @b_star.setter
    def b_star(self, value):
        self._direct_lattice.b = value

    @property
    def c_star(self):
        """Length of lattice vector :math:`\\mathbf{c^*}`."""
        return self._direct_lattice.c

    @c_star.setter
    def c_star(self, value):
        self._direct_lattice.c = value

    @property
    def alpha_star(self):
        """Angle between lattice vectors :math:`\\mathbf{b}^{\\ast}` and \
        :math:`\\mathbf{c}^{\\ast}` in **degrees**."""
        return self._direct_lattice.alpha

    @property
    def beta_star(self):
        """Angle between lattice vectors :math:`\\mathbf{c}^{\\ast}` and \
        :math:`\\mathbf{a}^{\\ast}` in **degrees**."""
        return self._direct_lattice.beta

    @property
    def gamma_star(self):
        """Angle between lattice vectors :math:`\\mathbf{a}^{\\ast}` and \
        :math:`\\mathbf{b}^{\\ast}` in **degrees**."""
        return self._direct_lattice.gamma

    @property
    def b1(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_1=\\mathbf{a}^*`."""
        return self._direct_lattice.a1

    @property
    def b2(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^*`."""
        return self._direct_lattice.a2

    @property
    def b3(self):
        """Reciprocal lattice vector :math:`\\mathbf{b}_3=\\mathbf{c}^*`."""
        return self._direct_lattice.a3

    @property
    def cos_alpha_star(self):
        """:math:`\\cos\\alpha^*`"""
        return self._direct_lattice.cos_alpha

    @property
    def cos_beta_star(self):
        """:math:`\\cos\\beta^*`"""
        return self._direct_lattice.cos_beta

    @property
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^*`"""
        return self._direct_lattice.cos_gamma

    @property
    def sin_alpha_star(self):
        """:math:`\\sin\\alpha^*`"""
        return self._direct_lattice.sin_alpha

    @property
    def sin_beta_star(self):
        """:math:`\\sin\\beta^*`"""
        return self._direct_lattice.sin_beta

    @property
    def sin_gamma_star(self):
        """:math:`\\sin\\gamma^*`"""
        return self._direct_lattice.sin_gamma

    @property
    def reciprocal_lattice(self):
        """Reciprocal lattice of this `Reciprocal3DLattice`."""
        return Crystal3DLattice(cell_matrix=np.linalg.inv(self.cell_matrix))

    @classmethod
    def triclinic(cls, a_star, b_star, c_star,
                  alpha_star, beta_star, gamma_star):
        """Generate a triclinic 3D lattice with lattice parameters \
            :math:`a^*, b^*, c^*, \\alpha^*, \\beta^*, \\gamma^*`."""
        return cls(a_star=a_star, b_star=b_star, c_star=c_star,
                   alpha_star=alpha_star, beta_star=beta_star,
                   gamma_star=gamma_star)

    @classmethod
    def monoclinic(cls, a_star, b_star, c_star, beta_star):
        """Generate a monoclinic 3D lattice with lattice parameters \
            :math:`a^*, b^*, c^*, \\beta^*`."""
        return cls(a_star=a_star, b_star=b_star, c_star=c_star,
                   alpha_star=90, beta_star=beta_star, gamma_star=90)

    @classmethod
    def orthorhombic(cls, a_star, b_star, c_star):
        """Generate an orthorhombic 3D lattice with lattice parameters \
            :math:`a^*, b^*, c^*`."""
        return cls(a_star=a_star, b_star=b_star, c_star=c_star,
                   alpha_star=90, beta_star=90, gamma_star=90)

    @classmethod
    def tetragonal(cls, a_star, c_star):
        """Generate a tetragonal 3D lattice with lattice parameters \
            :math:`a^*, c^*`."""
        return cls(a_star=a_star, b_star=a_star, c_star=c_star,
                   alpha_star=90, beta_star=90, gamma_star=90)

    @classmethod
    def hexagonal(cls, a_star, c_star):
        """Generate a hexagonal 3D lattice with lattice parameters \
            :math:`a^*, c^*`."""
        return cls(a_star=a_star, b_star=a_star, c_star=c_star,
                   alpha_star=90, beta_star=90, gamma_star=120)

    @classmethod
    def cubic(cls, a_star):
        """Generate a cubic 3D lattice with lattice parameter \
            :math:`a^*`."""
        return cls(a_star=a_star, b_star=a_star, c_star=a_star,
                   alpha_star=90, beta_star=90, gamma_star=90)

ReciprocalLattice = Reciprocal3DLattice


# class Bravais3DLattice:
#     """Class for bravais lattices."""
#     def __init__(self, crystal_system=None, centering=None,
#                  symbol=None):
#         pass

# BravaisLattice = Bravais3DLattice


# class SimpleCubicLattice(BravaisLattice):
#     """Abstract representation of simple cubic lattice."""
#     lattice_points = [[0.0, 0.0, 0.0]]


# class BodyCenteredCubicLattice(BravaisLattice):
#     """Abstract representation of body-centered cubic lattice."""
#     lattice_points = [[0.0, 0.0, 0.0],
#                       [0.5, 0.5, 0.5]]


# class FaceCenteredCubicLattice(BravaisLattice):
#     """Abstract representation of face-centered cubic lattice."""
#     lattice_points = [[0.0, 0.0, 0.0],
#                       [0.5, 0.5, 0.0],
#                       [0.0, 0.5, 0.5],
#                       [0.5, 0.0, 0.5]]
