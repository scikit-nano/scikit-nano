# -*- coding: utf-8 -*-
"""
=============================================================================
Crystal lattice classes (:mod:`sknano.core.crystallography._xtal_lattices`)
=============================================================================

.. currentmodule:: sknano.core.crystallography._xtal_lattices

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty
from functools import total_ordering

import numpy as np

from sknano.core import BaseClass
from sknano.core.math import Vector, Point, zhat, rotation_matrix

__all__ = ['LatticeBase', 'ReciprocalLatticeBase',
           'DirectLatticeMixin', 'ReciprocalLatticeMixin',
           'Direct2DLatticeMixin', 'Direct3DLatticeMixin',
           'Reciprocal2DLatticeMixin', 'Reciprocal3DLatticeMixin',
           'Crystal2DLattice', 'Reciprocal2DLattice',
           'CrystalLattice', 'ReciprocalLattice',
           'Crystal3DLattice', 'Reciprocal3DLattice']
# __all__ += ['BravaisLattice', 'Bravais3DLattice',
#            'SimpleCubicLattice', 'BodyCenteredCubicLattice',
#            'FaceCenteredCubicLattice']


@total_ordering
class LatticeBase(BaseClass):
    """Base class for crystallographic lattice objects.

    Parameters
    ----------
    nd : int
    cell_matrix : array_like
    orientation_matrix : array_like, optional

    """

    def __init__(self, nd=None, cell_matrix=None, orientation_matrix=None):
        super().__init__()

        self.nd = nd
        self.offset = Point()
        if cell_matrix is not None and orientation_matrix is None:
            orientation_matrix = cell_matrix.T * self.fractional_matrix

        if orientation_matrix is None:
            orientation_matrix = np.asmatrix(np.identity(3))

        self.orientation_matrix = np.asmatrix(orientation_matrix)
        self.lattice_type = None

    def __dir__(self):
        return ['nd', 'offset', 'orientation_matrix']

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self is other or \
                all([np.allclose(getattr(self, attr), getattr(other, attr))
                     for attr in dir(self)])

    def __lt__(self, other):
        if isinstance(other, type(self)):
            try:
                return self.cell_volume < other.cell_volume
            except AttributeError:
                return self.cell_area < other.cell_area

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
            :attr:`~sknano.core.crystallography.LatticeBase.cell_matrix`."""
        return self.cell_matrix

    @property
    def fractional_matrix(self):
        """Transformation matrix to convert from cartesian coordinates to \
            fractional coordinates."""
        return np.linalg.inv(self.ortho_matrix)

    @property
    def metric_tensor(self):
        """Metric tensor."""
        return self.cell_matrix * self.cell_matrix.T

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

    def wrap_fractional_coordinate(self, p, epsilon=1e-6, pbc=None):
        """Wrap fractional coordinate to lie within unit cell.

        Parameters
        ----------
        p : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        if pbc is None:
            pbc = np.asarray(np.ones(3), dtype=bool)

        p = np.ma.array(p, mask=~pbc)
        p = np.ma.fmod(p, 1)
        p[np.ma.where(p < 0)] += 1
        p[np.ma.where(p > 1 - epsilon)] -= 1
        p[np.ma.where(np.logical_or((p > 1 - epsilon), (p < epsilon)))] = 0
        p.mask = np.ma.nomask
        return p.tolist()

    def wrap_cartesian_coordinate(self, p, pbc=None):
        """Wrap cartesian coordinate to lie within unit cell.

        Parameters
        ----------
        p : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        return self.fractional_to_cartesian(
            self.wrap_fractional_coordinate(self.cartesian_to_fractional(p),
                                            pbc=pbc))

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
                np.asmatrix(
                    rotation_matrix(angle=angle, axis=axis,
                                    anchor_point=anchor_point,
                                    rot_point=rot_point,
                                    from_vector=from_vector,
                                    to_vector=to_vector, degrees=degrees,
                                    verbose=verbose, **kwargs))
            # print('transform_matrix: {}'.format(transform_matrix))

            # transform_matrix = \
            #     transformation_matrix(angle=angle, axis=axis,
            #                           anchor_point=anchor_point,
            #                           rot_point=rot_point,
            #                           from_vector=from_vector,
            #                           to_vector=to_vector, degrees=degrees,
            #                           verbose=verbose, **kwargs)

        self.orientation_matrix = \
            transform_matrix * self.orientation_matrix

    def translate(self, t):
        """Translate lattice.

        Parameters
        ----------
        t : :class:`Vector`

        See Also
        --------
        core.math.translate

        """
        self.offset.translate(t)


class ReciprocalLatticeBase(LatticeBase):
    """Base class for crystallographic reciprocal lattice objects.

    Parameters
    ----------
    direct_lattice : :class:`Crystal2DLattice` or :class:`Crystal3DLattice`
    nd : int
    """
    def __init__(self, direct_lattice, nd):
        self._direct_lattice = direct_lattice
        super().__init__(
            nd=nd, cell_matrix=self._direct_lattice.cell_matrix,
            orientation_matrix=self._direct_lattice.orientation_matrix)

    def __getattr__(self, name):
        if name != '_direct_lattice':
            return getattr(self._direct_lattice, name)


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
            (self.cos_gamma_star * self.cos_alpha_star - self.cos_beta_star) /
            (self.sin_gamma_star * self.sin_alpha_star), decimals=10)

    @property
    def cos_gamma(self):
        """:math:`\\cos\\gamma`"""
        print('calling Direct3DLatticeMixin.cos_gamma')
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


class Crystal2DLattice(LatticeBase, Reciprocal2DLatticeMixin):
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
        self._a = float(value)
        self._update_ortho_matrix()

    @property
    def b(self):
        """Length of lattice_vector :math:`\\mathbf{b}`."""
        return np.around(self._b, decimals=10)

    @b.setter
    def b(self, value):
        self._b = float(value)
        self._update_ortho_matrix()

    @property
    def gamma(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{b}` in **degrees**."""
        try:
            return np.around(float(self._gamma), decimals=10)
        except TypeError:
            return None

    @gamma.setter
    def gamma(self, value):
        self._gamma = float(value)
        self._update_ortho_matrix()

    @property
    def a1(self):
        """Lattice vector :math:`\\mathbf{a}_1=\\mathbf{a}`."""
        return Vector(self.cell_matrix[0, :].A.flatten())[:2]

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
        """Crystal cell area."""
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
        self._ortho_matrix = np.matrix([[m11, m12, 0.0],
                                        [0.0, m22, 0.0],
                                        [0.0, 0.0, 1.0]])

    @property
    def reciprocal_lattice(self):
        """Return `Crystal2DLattice` reciprocal lattice."""
        return Reciprocal2DLattice(
            cell_matrix=np.linalg.inv(self.cell_matrix).T)

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

    @property
    def lengths(self):
        """Tuple of lattice parameter lengths :math:`a, b`."""
        return self.a, self.b

    @property
    def angles(self):
        """Lattice parameter angle \\gamma`."""
        return self.gamma

    @property
    def lattice_parameters(self):
        """Tuple of lattice parameters `a`, `b`, `gamma`."""
        return self.a, self.b, self.gamma

    @property
    def lengths_and_angles(self):
        """Alias for attr:`Crystal2DLattice.lattice_parameters`."""
        return self.lattice_parameters

    @property
    def lattice_vectors(self):
        """Tuple of lattice vectors :math:`\\mathbf{a}_1, \\mathbf{a}_2`."""
        return self.a1, self.a2


class Reciprocal2DLattice(ReciprocalLatticeBase, Direct2DLatticeMixin):
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

    @gamma_star.setter
    def gamma_star(self, value):
        self._direct_lattice.gamma = value

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
        return Crystal2DLattice(cell_matrix=np.linalg.inv(self.cell_matrix).T)

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

    @property
    def lengths(self):
        """Tuple of lattice parameter lengths :math:`a^*, b^*`."""
        return self.a_star, self.b_star

    @property
    def angles(self):
        """Lattice parameter angle \\gamma^*`."""
        return self.gamma_star

    @property
    def lattice_parameters(self):
        """Tuple of lattice parameters `a^*`, `b^*`, `gamma^*`."""
        return self.a_star, self.b_star, self.gamma_star

    @property
    def lattice_vectors(self):
        """Tuple of lattice vectors :math:`\\mathbf{b}_1, \\mathbf{b}_2`."""
        return self.b1, self.b2


class Crystal3DLattice(LatticeBase, ReciprocalLatticeMixin):
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

        if all([v is not None for v in (a1, a2, a3)]):
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

        if all([p is not None for p in (a, b, c, alpha, beta, gamma)]):
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
        self._a = float(value)
        self._update_ortho_matrix()

    @property
    def b(self):
        """Length of lattice_vector :math:`\\mathbf{b}`."""
        return np.around(self._b, decimals=10)

    @b.setter
    def b(self, value):
        self._b = float(value)
        self._update_ortho_matrix()

    @property
    def c(self):
        """Length of lattice vector :math:`\\mathbf{c}`."""
        return np.around(self._c, decimals=10)

    @c.setter
    def c(self, value):
        self._c = float(value)
        self._update_ortho_matrix()

    @property
    def alpha(self):
        """Angle between lattice vectors :math:`\\mathbf{b}` and \
        :math:`\\mathbf{c}` in **degrees**."""
        return np.around(self._alpha, decimals=10)

    @alpha.setter
    def alpha(self, value):
        self._alpha = float(value)
        self._update_ortho_matrix()

    @property
    def beta(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{c}` in **degrees**."""
        return np.around(self._beta, decimals=10)

    @beta.setter
    def beta(self, value):
        self._beta = float(value)
        self._update_ortho_matrix()

    @property
    def gamma(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{b}` in **degrees**."""
        return np.around(self._gamma, decimals=10)

    @gamma.setter
    def gamma(self, value):
        self._gamma = float(value)
        self._update_ortho_matrix()

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
        """Crystal cell volume."""
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
        return Reciprocal3DLattice(
            cell_matrix=np.linalg.inv(self.cell_matrix).T)

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

    @property
    def lengths(self):
        """Tuple of lattice vector lengths :math:`a, b, c`."""
        return self.a, self.b, self.c

    @property
    def angles(self):
        """Tuple of lattice parameter angles \
            :math:`\\alpha, \\beta, \\gamma`."""
        return self.alpha, self.beta, self.gamma

    @property
    def lattice_parameters(self):
        """Tuple of lattice parameters \
            `a`, `b`, `c`, `alpha`, `beta`, `gamma`."""
        return self.a, self.b, self.c, self.alpha, self.beta, self.gamma

    @property
    def lengths_and_angles(self):
        """Alias for attr:`Crystal3DLattice.lattice_parameters`."""
        return self.lattice_parameters

    @property
    def lattice_vectors(self):
        """Tuple of lattice vectors \
            :math:`\\mathbf{a}_1, \\mathbf{a}_2, \\mathbf{a}_3`."""
        return self.a1, self.a2, self.a3


CrystalLattice = Crystal3DLattice


class Reciprocal3DLattice(ReciprocalLatticeBase, DirectLatticeMixin):
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

    @alpha_star.setter
    def alpha_star(self, value):
        self._direct_lattice.alpha = value

    @property
    def beta_star(self):
        """Angle between lattice vectors :math:`\\mathbf{c}^{\\ast}` and \
            :math:`\\mathbf{a}^{\\ast}` in **degrees**."""
        return self._direct_lattice.beta

    @beta_star.setter
    def beta_star(self, value):
        self._direct_lattice.beta = value

    @property
    def gamma_star(self):
        """Angle between lattice vectors :math:`\\mathbf{a}^{\\ast}` and \
            :math:`\\mathbf{b}^{\\ast}` in **degrees**."""
        return self._direct_lattice.gamma

    @gamma_star.setter
    def gamma_star(self, value):
        self._direct_lattice.gamma = value

    @property
    def b1(self):
        """Lattice vector :math:`\\mathbf{b}_1=\\mathbf{a}^{\\ast}`."""
        return self._direct_lattice.a1

    @property
    def b2(self):
        """Lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^{\\ast}`."""
        return self._direct_lattice.a2

    @property
    def b3(self):
        """Lattice vector :math:`\\mathbf{b}_3=\\mathbf{c}^{\\ast}`."""
        return self._direct_lattice.a3

    @property
    def cos_alpha_star(self):
        """:math:`\\cos\\alpha^{\\ast}`"""
        return self._direct_lattice.cos_alpha

    @property
    def cos_beta_star(self):
        """:math:`\\cos\\beta^{\\ast}`"""
        return self._direct_lattice.cos_beta

    @property
    def cos_gamma_star(self):
        """:math:`\\cos\\gamma^{\\ast}`"""
        return self._direct_lattice.cos_gamma

    @property
    def sin_alpha_star(self):
        """:math:`\\sin\\alpha^{\\ast}`"""
        return self._direct_lattice.sin_alpha

    @property
    def sin_beta_star(self):
        """:math:`\\sin\\beta^{\\ast}`"""
        return self._direct_lattice.sin_beta

    @property
    def sin_gamma_star(self):
        """:math:`\\sin\\gamma^{\\ast}`"""
        return self._direct_lattice.sin_gamma

    @property
    def reciprocal_lattice(self):
        """Reciprocal lattice of this `Reciprocal3DLattice`."""
        return Crystal3DLattice(cell_matrix=np.linalg.inv(self.cell_matrix).T)

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

    def todict(self):
        """Return `dict` of `Reciprocal3DLattice` parameters."""
        return dict(a_star=self.a_star, b_star=self.b_star,
                    c_star=self.c_star, alpha_star=self.alpha_star,
                    beta_star=self.beta_star, gamma_star=self.gamma_star)

    @property
    def lengths(self):
        """Tuple of lattice vector lengths :math:`a^*, b^*, c^*`."""
        return self.a_star, self.b_star, self.c_star

    @property
    def angles(self):
        """Tuple of lattice parameter angles \
            :math:`\\alpha^*, \\beta^*, \\gamma^*`."""
        return self.alpha_star, self.beta_star, self.gamma_star

    @property
    def lattice_parameters(self):
        """Tuple of lattice parameters \
            `a^*`, `b^*`, `c^*`, `alpha^*`, `beta^*`, `gamma^*`."""
        return self.a_star, self.b_star, self.c_star, \
            self.alpha_star, self.beta_star, self.gamma_star

    @property
    def lattice_vectors(self):
        """Tuple of lattice vectors \
            :math:`\\mathbf{b}_1, \\mathbf{b}_2, \\mathbf{b}_3`."""
        return self.b1, self.b2, self.b3


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
