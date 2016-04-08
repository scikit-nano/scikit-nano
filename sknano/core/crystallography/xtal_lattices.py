# -*- coding: utf-8 -*-
"""
=============================================================================
Crystal lattice classes (:mod:`sknano.core.crystallography.xtal_lattices`)
=============================================================================

.. currentmodule:: sknano.core.crystallography.xtal_lattices

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty
from functools import total_ordering

import numpy as np

# from tabulate import tabulate

from sknano.core import BaseClass, TabulateMixin
from sknano.core.geometric_regions import Cuboid, Parallelepiped, \
    generate_bounding_box
from sknano.core.math import Vector, Point, zhat, transformation_matrix
from .extras import pbc_diff

__all__ = ['LatticeBase', 'ReciprocalLatticeBase',
           'DirectLatticeMixin', 'ReciprocalLatticeMixin',
           'Direct2DLatticeMixin', 'Direct3DLatticeMixin',
           'Reciprocal2DLatticeMixin', 'Reciprocal3DLatticeMixin',
           'Crystal2DLattice', 'Reciprocal2DLattice',
           'Crystal3DLattice', 'Reciprocal3DLattice',
           'Domain', 'generate_lattice']
# __all__ += ['BravaisLattice', 'Bravais3DLattice',
#            'SimpleCubicLattice', 'BodyCenteredCubicLattice',
#            'FaceCenteredCubicLattice']


class Domain(TabulateMixin, BaseClass):
    """Container class for molecular dynamics simulation box metadata."""
    def __init__(self):
        """
        Attributes
        ----------
        bounding_box : :class:`~sknano.core.geometric_regions.Cuboid`
        triclinic : :class:`~python:bool`
        xy, xz, yz : :class:`~python:float`
        lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`.

        """
        self.bounding_box = Cuboid()
        self.triclinic = False
        self.xy = self.xz = self.yz = 0.0
        self.lattice = None

    def __str__(self):
        strrep = self._table_title_str()
        attrs = ['triclinic', 'xy', 'xz', 'yz']
        values = [getattr(self, attr) for attr in attrs]
        table = self._tabulate(list(zip(attrs, values)))
        strrep = '\n'.join((strrep, table))
        lattice = self.lattice
        if lattice is not None:
            strrep = '\n'.join((strrep, str(lattice)))
        return strrep

    @property
    def tilt_factors(self):
        """Domain tilt factors `xy`, `xz`, `yz`."""
        return self.xy, self.xz, self.yz

    @property
    def xlo(self):
        """Alias for :attr:`Domain.bounding_box.xmin`"""
        return self.bounding_box.xmin

    @xlo.setter
    def xlo(self, value):
        self.bounding_box.xmin = value

    @property
    def xhi(self):
        """Alias for :attr:`Domain.bounding_box.xmax`"""
        return self.bounding_box.xmax

    @xhi.setter
    def xhi(self, value):
        self.bounding_box.xmax = value

    @property
    def ylo(self):
        """Alias for :attr:`Domain.bounding_box.ymin`"""
        return self.bounding_box.ymin

    @ylo.setter
    def ylo(self, value):
        self.bounding_box.ymin = value

    @property
    def yhi(self):
        """Alias for :attr:`Domain.bounding_box.ymax`"""
        return self.bounding_box.ymax

    @yhi.setter
    def yhi(self, value):
        self.bounding_box.ymax = value

    @property
    def zlo(self):
        """Alias for :attr:`Domain.bounding_box.zmin`"""
        return self.bounding_box.zmin

    @zlo.setter
    def zlo(self, value):
        self.bounding_box.zmin = value

    @property
    def zhi(self):
        """Alias for :attr:`Domain.bounding_box.zmax`"""
        return self.bounding_box.zmax

    @zhi.setter
    def zhi(self, value):
        self.bounding_box.zmax = value

    @property
    def lengths(self):
        """:class:`~python:tuple` of side lengths"""
        return self.bounding_box.lengths

    @property
    def lx(self):
        """Alias for :attr:`Domain.bounding_box.lx`."""
        return self.bounding_box.lx

    @property
    def ly(self):
        """Alias for :attr:`Domain.bounding_box.ly`."""
        return self.bounding_box.ly

    @property
    def lz(self):
        """Alias for :attr:`Domain.bounding_box.lz`."""
        return self.bounding_box.lz

    @property
    def xlo_bound(self):
        """Triclinic bounding box minimum extent in the x-dimension"""
        xlo = self.xlo
        xy = self.xy
        xz = self.xz
        return xlo + min((0.0, xy, xz, xy + xz))

    @xlo_bound.setter
    def xlo_bound(self, value):
        xy = self.xy
        xz = self.xz
        self.xlo = value - min((0.0, xy, xz, xy + xz))

    @property
    def xhi_bound(self):
        """Triclinic bounding box maximum extent in the x-dimension"""
        xhi = self.xhi
        xy = self.xy
        xz = self.xz
        return xhi + max((0.0, xy, xz, xy + xz))

    @xhi_bound.setter
    def xhi_bound(self, value):
        xy = self.xy
        xz = self.xz
        self.xhi = value - max((0.0, xy, xz, xy + xz))

    @property
    def ylo_bound(self):
        """Triclinic bounding box minimum extent in the y-dimension"""
        return self.ylo + min((0.0, self.yz))

    @ylo_bound.setter
    def ylo_bound(self, value):
        self.ylo = value - min((0.0, self.yz))

    @property
    def yhi_bound(self):
        """Triclinic bounding box maximum extent in the y-dimension"""
        return self.yhi + max((0.0, self.yz))

    @yhi_bound.setter
    def yhi_bound(self, value):
        self.yhi = value - max((0.0, self.yz))

    @property
    def zlo_bound(self):
        """Triclinic bounding box minimum extent in the z-dimension"""
        return self.zlo

    @zlo_bound.setter
    def zlo_bound(self, value):
        self.zlo = value

    @property
    def zhi_bound(self):
        """Triclinic bounding box maximum extent in the z-dimension"""
        return self.zhi

    @zhi_bound.setter
    def zhi_bound(self, value):
        self.zhi = value

    def update(self, from_lattice=None, from_region=None, from_array=None,
               allow_triclinic_box=False, pad_box=False,
               pad_tol=0.01, xpad=10., ypad=10., zpad=10., verbose=False,
               **kwargs):
        """Update simulation domain attributes from lattice."""
        bounding_box = \
            generate_bounding_box(from_lattice=from_lattice,
                                  from_region=from_region,
                                  from_array=from_array,
                                  verbose=verbose)

        if pad_box and from_array is not None:
            coords = from_array
            boxpad = {'x': xpad, 'y': ypad, 'z': zpad}
            # for dim, pad in boxpad.items():
            for i, dim in enumerate(('x', 'y', 'z')):
                pad = boxpad[dim]
                dmin = dim + 'min'
                dmax = dim + 'max'
                if abs(getattr(bounding_box, dmin) -
                       coords[:, i].min()) < pad - pad_tol:
                    setattr(bounding_box, dmin,
                            getattr(bounding_box, dmin) - pad)
                if abs(getattr(bounding_box, dmax) -
                       coords[:, i].max()) < pad - pad_tol:
                    setattr(bounding_box, dmax,
                            getattr(bounding_box, dmax) + pad)

        if allow_triclinic_box and from_lattice is not None:
            self.lattice = lattice = from_lattice
            if not np.allclose(np.radians(lattice.angles),
                               np.pi / 2 * np.ones(3)):
                self.triclinic = True
                a, b, c = lattice.lengths
                cos_alpha, cos_beta, cos_gamma = \
                    np.cos(np.radians(lattice.angles))
                self.xy = xy = b * cos_gamma
                self.xz = xz = c * cos_beta
                self.yz = \
                    (b * c * cos_alpha - xy * xz) / np.sqrt(b ** 2 - xy ** 2)

                self.xlo_bound = bounding_box.xmin
                self.xhi_bound = bounding_box.xmax
                self.ylo_bound = bounding_box.ymin
                self.yhi_bound = bounding_box.ymax
                self.zlo_bound = bounding_box.zmin
                self.zhi_bound = bounding_box.zmax
        else:
            self.bounding_box = bounding_box

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict()


def generate_lattice(from_domain=None, from_region=None, offset=None):
    """Return a :class:`~sknano.core.crystallography.Crystal3DLattice`.

    Parameters
    ----------
    from_domain : :class:`~sknano.core.crystallography.Domain`
    from_region : :class:`~sknano.core.geometric_regions.Geometric3DRegion`
    verbose : :class:`~python:bool`

    Returns
    -------
    lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`

    """
    if all([obj is None for obj in (from_domain, from_region)]):
        return None

    if from_domain is not None:
        domain = from_domain
        lx, ly, lz = domain.lengths
        xy, xz, yz = domain.tilt_factors

        a = lx
        b = np.sqrt(ly ** 2 + xy ** 2)
        c = np.sqrt(lz ** 2 + xz ** 2 + yz ** 2)
        alpha = np.degrees(np.arccos((xy * xz + ly * yz) / (b * c)))
        beta = np.degrees(np.arccos(xz / c))
        gamma = np.degrees(np.arccos(xy / b))

        lattice = \
            Crystal3DLattice(a=a, b=b, c=c, alpha=alpha, beta=beta,
                             gamma=gamma, offset=offset)
        return lattice


@total_ordering
class LatticeBase(TabulateMixin, BaseClass):
    """Base class for crystallographic lattice objects.

    Parameters
    ----------
    nd : int
    cell_matrix : array_like
    orientation_matrix : array_like, optional
    offset : array_like, optional

    """
    def __init__(self, nd=None, cell_matrix=None, orientation_matrix=None,
                 offset=None):
        super().__init__()

        self.nd = nd
        if cell_matrix is not None and orientation_matrix is None:
            orientation_matrix = cell_matrix.T * self.fractional_matrix

        if orientation_matrix is None:
            orientation_matrix = np.asmatrix(np.identity(3))

        self.orientation_matrix = np.asmatrix(orientation_matrix)
        self.lattice_type = None
        self._offset = Point(offset, nd=3)
        self._update_region()

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
    def offset(self):
        """Lattice offset."""
        return self._offset

    @offset.setter
    def offset(self, value):
        # self.translate(Vector(p0=self._offset, p=Point(value)))
        self._offset[:] = Point(value)

    @property
    def ortho_matrix(self):
        """Transformation matrix to convert from fractional coordinates to \
            cartesian coordinates.

        .. math::

           [M] = \\begin{pmatrix}
           \\mathbf{a} & \\mathbf{b} & \\mathbf{c}
           \\end{pmatrix}

           =\\begin{pmatrix}
           a & b\\cos\\gamma & c\\cos\\beta\\\\
           0 & b\\sin\\gamma & c(\\cos\\alpha - \\cos\\beta\\cos\\gamma)/\\sin\\gamma\\\\
           0 & 0 & c\\sin\\alpha\\sin\\beta\\sin\\gamma^*/\\sin\\gamma
           \\end{pmatrix}

        """
        return self._ortho_matrix

    @property
    def cell_matrix(self):
        """Matrix of lattice row vectors.

        Same as the matrix transpose of the matrix product
        :attr:`LatticeBase.orientation_matrix` :math:`\\times`
        :attr:`LatticeBase.ortho_matrix`:

        .. math::

           [A] = ([R][M])^T
           \\begin{pmatrix}
           \\mathbf{a}\\\\
           \\mathbf{b}\\\\
           \\mathbf{c}
           \\end{pmatrix}
           =([R][\\mathbf{a}\\,\\mathbf{b}\\,\\mathbf{c}])^T

        """
        return (self.orientation_matrix * self.ortho_matrix).T

    @property
    def cell(self):
        """Alias for \
            :attr:`~sknano.core.crystallography.LatticeBase.cell_matrix`."""
        return self.cell_matrix

    @property
    def matrix(self):
        """Alias for \
            :attr:`~sknano.core.crystallography.LatticeBase.cell_matrix`."""
        return self.cell_matrix

    @property
    def fractional_matrix(self):
        """Transformation matrix to convert from cartesian coordinates to \
            fractional coordinates.

        The fractional matrix :math:`[Q]` is given by the inverse
        of the orthogonal matrix :math:`[M]`
        (see: :attr:`~LatticeBase.ortho_matrix`),
        i.e. :math:`[Q]=[M]^{-1}`

        .. math::

           [Q] = [M]^{-1} =\\begin{pmatrix}
           1/a & -1/a\\tan\\gamma & (\\cos\\alpha\\cos\\gamma - \\cos\\beta)/
           a\\phi\\sin\\gamma\\\\
           0 & 1/b\\sin\\gamma & (\\cos\\beta\\cos\\gamma - \\cos\\alpha)/
           b\\phi\\sin\\gamma\\\\
           0 & 0 & \\sin\\gamma/c\\phi
           \\end{pmatrix}

        where:

        .. math::

           \\begin{align*}
           \\phi &= \\frac{V}{abc}\\\\
                 &=\\sqrt{1 - \\cos^2\\alpha - \\cos^2\\beta - \\cos^2\\gamma
                           + 2\\cos\\alpha\\cos\\beta\\cos\\gamma}\\\\
                 &=\\sin\\alpha\\sin\\beta\\sin\\gamma^*
           \\end{align*}

        where :math:`V` is the volume of the unit cell.

        """
        return np.linalg.inv(self.ortho_matrix)

    @property
    def metric_tensor(self):
        """Metric tensor."""
        return self.cell_matrix * self.cell_matrix.T

    @property
    def bounding_box(self):
        """:attr:`~LatticeBase.region` :attr:`~Parallelepiped.bounding_box`."""
        return self.region.bounding_box

    @property
    def region(self):
        """:class:`Parallelepiped` defined by lattice vectors."""
        try:
            return self._region
        except AttributeError:
            self._update_region()
            return self._region

    def _update_region(self):
        cell_matrix = self.cell_matrix
        o = self.offset
        u, v, w = \
            map(Vector, [cell_matrix[ri].A.flatten() for ri in range(3)])
        # u, v, w = \
        #     map(Vector, [ortho_matrix[:, ri].A.flatten() for ri in range(3)])
        self._region = Parallelepiped(o, u, v, w)

    @property
    def centroid(self):
        """Region centroid."""
        return self.region.centroid

    def fractional_diff(self, fcoords1, fcoords2):
        """Compute difference between fractional coordinates.

        See Also
        --------
        core.crystallography.pbc_diff

        """
        return pbc_diff(fcoords1, fcoords2)

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

    def wrap_fractional_coordinate(self, p, epsilon=1e-8, pbc=None):
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
        return np.asarray(p.tolist())

    def wrap_fractional_coordinates(self, points, epsilon=1e-8, pbc=None):
        """Wrap array of fractional coordinates to lie within unit cell.

        Parameters
        ----------
        p : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        return np.asarray([self.wrap_fractional_coordinate(p, epsilon=epsilon,
                                                           pbc=pbc)
                           for p in points])

    def wrap_cartesian_coordinate(self, p, epsilon=1e-8, pbc=None):
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
                                            epsilon=epsilon, pbc=pbc))

    def wrap_cartesian_coordinates(self, points, epsilon=1e-8, pbc=None):
        """Wrap array of cartesian coordinates to lie within unit cell.

        Parameters
        ----------
        p : array_like

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        return np.asarray([self.wrap_cartesian_coordinate(p, epsilon=epsilon,
                                                          pbc=pbc)
                           for p in points])

    def rotate(self, **kwargs):
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
        if self.nd == 2 and kwargs.get('axis', None) is None:
            kwargs['axis'] = 'z'
        if kwargs.get('anchor_point', None) is None:
            kwargs['anchor_point'] = self.offset

        if kwargs.get('transform_matrix', None) is None:
            kwargs['transform_matrix'] = transformation_matrix(**kwargs)
        transform_matrix = kwargs['transform_matrix']
        if transform_matrix.shape == (4, 4):
            # tvec = translation_from_matrix(transform_matrix)
            # print('translation_part: {}'.format(tvec))
            # self.translate(tvec)
            transform_matrix = transform_matrix[:3, :3]
        transform_matrix = np.asmatrix(transform_matrix)
        self.orientation_matrix = transform_matrix * self.orientation_matrix

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
        self.region.translate(t)


class ReciprocalLatticeBase(LatticeBase):
    """Base class for crystallographic reciprocal lattice objects.

    Parameters
    ----------
    direct_lattice : :class:`Crystal2DLattice` or :class:`Crystal3DLattice`
    nd : int
    """
    def __init__(self, direct_lattice, nd, offset=None):
        self._direct_lattice = direct_lattice
        super().__init__(
            nd=nd, cell_matrix=self._direct_lattice.cell_matrix,
            orientation_matrix=self._direct_lattice.orientation_matrix,
            offset=offset)

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
        # b2.p0 = self.offset
        return b2.cross(zhat) / self.cell_area

    @property
    def a2(self):
        """2D lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        b1 = Vector()
        b1[:2] = self.b1
        # z = zhat.copy()
        # z.p0 = self.offset
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
        # a2.p0 = self.offset
        return a2.cross(zhat)[:2] / self.cell_area

    @property
    def b2(self):
        """2D reciprocal lattice vector :math:`\\mathbf{b}_2=\\mathbf{b}^{*}`.
        """
        a1 = Vector()
        a1[:2] = self.a1
        # z = zhat.copy()
        # z.p0 = self.offset
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

    def __init__(self, a=None, b=None, gamma=None, a1=None, a2=None,
                 cell_matrix=None, orientation_matrix=None, offset=None):

        if cell_matrix is not None:
            cell_matrix = np.asarray(cell_matrix)
            a1 = np.array(cell_matrix[0])
            a2 = np.array(cell_matrix[1])

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
                         orientation_matrix=orientation_matrix,
                         offset=offset)

        self.fmtstr = "a={a!r}, b={b!r}, gamma={gamma!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['a', 'b', 'gamma'])
        return attrs

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        lattice_params = \
            self._tabulate(list(zip(('a', 'b', 'γ'),
                                    self.lengths_and_angles)))
        title = '.'.join((objstr, 'lengths_and_angles'))
        strrep = '\n'.join((strrep, title, lattice_params))
        cosines = np.around(np.cos(np.radians(self.angles)), decimals=2)
        cosines = self._tabulate(list(zip(('cos(γ)',), (cosines,))))
        title = '.'.join((objstr, 'offset'))
        offset = self._tabulate([['Lattice offset', self.offset]])
        strrep = '\n'.join((strrep, title, offset))

        lattice_vectors = \
            self._tabulate(list(zip(('a1', 'a2'), self.lattice_vectors)))
        title = '.'.join((objstr, 'lattice_vectors'))
        strrep = '\n'.join((strrep, title, lattice_vectors))

        region = self.region
        title = '.'.join((objstr, region.__class__.__qualname__))
        strrep = '\n'.join((strrep, title, str(region)))

        bbox = region.bounding_box
        title = '.'.join((title, bbox.__class__.__qualname__))
        strrep = '\n'.join((strrep, title, str(bbox)))

        return strrep

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
        # return Vector(self.ortho_matrix[:, 0].A.flatten(),
        #               p0=self.offset)[:2]
        return Vector(self.cell_matrix[0].A.flatten(), p0=self.offset)[:2]

    @property
    def a2(self):
        """Lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        # return Vector(self.ortho_matrix[:, 1].A.flatten(),
        #               p0=self.offset)[:2]
        return Vector(self.cell_matrix[1].A.flatten(), p0=self.offset)[:2]

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
                 b1=None, b2=None, cell_matrix=None, orientation_matrix=None,
                 offset=None):

        direct_lattice = \
            Crystal2DLattice(a=a_star, b=b_star, gamma=gamma_star,
                             a1=b1, a2=b2, cell_matrix=cell_matrix,
                             orientation_matrix=orientation_matrix)
        super().__init__(direct_lattice, nd=2, offset=offset)

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
                 orientation_matrix=None, offset=None):

        if all([p is None for p in (a, b, c, alpha, beta, gamma,
                                    a1, a2, a3, cell_matrix)]):
            errmsg = 'Expected lattice parameters ' + \
                '`a`, `b`, `c`, `alpha`, `beta`, `gamma`\n' + \
                'or lattice vectors `a1`, `a2`, `a3`\n' + \
                'or a 3x3 matrix with rows of lattice vectors.'
            raise ValueError(errmsg)

        if cell_matrix is not None:
            cell_matrix = np.asarray(cell_matrix)
            a1 = np.array(cell_matrix[0])
            a2 = np.array(cell_matrix[1])
            a3 = np.array(cell_matrix[2])

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
                         orientation_matrix=orientation_matrix,
                         offset=offset)

        self.fmtstr = "a={a!r}, b={b!r}, c={c!r}, " + \
            "alpha={alpha!r}, beta={beta!r}, gamma={gamma!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['a', 'b', 'c', 'alpha', 'beta', 'gamma'])
        return attrs

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        lattice_params = \
            self._tabulate(list(zip(('a', 'b', 'c', 'α', 'β', 'γ'),
                                    self.lengths_and_angles)))
        title = '.'.join((objstr, 'lengths_and_angles'))
        strrep = '\n'.join((strrep, title, lattice_params))
        cosines = np.around(np.cos(np.radians(self.angles)), decimals=2)
        cosines = self._tabulate(list(zip(('cos(α)', 'cos(β)', 'cos(γ)'),
                                          cosines)))
        title = '.'.join((objstr, 'offset'))
        offset = self._tabulate([['Lattice offset', self.offset]])
        strrep = '\n'.join((strrep, title, offset))

        lattice_vectors = \
            self._tabulate(list(zip(('a1', 'a2', 'a3'),
                                    self.lattice_vectors)))
        title = '.'.join((objstr, 'lattice_vectors'))
        strrep = '\n'.join((strrep, title, lattice_vectors))

        region = self.region
        title = '.'.join((objstr, region.__class__.__qualname__))
        strrep = '\n'.join((strrep, title, str(region)))

        bbox = region.bounding_box
        title = '.'.join((title, bbox.__class__.__qualname__))
        strrep = '\n'.join((strrep, title, str(bbox)))

        return strrep

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
        # return Vector(self.ortho_matrix[:, 0].A.flatten(), p0=self.offset)
        return Vector(self.cell_matrix[0].A.flatten(), p0=self.offset)

    @property
    def a2(self):
        """Lattice vector :math:`\\mathbf{a}_2=\\mathbf{b}`."""
        # return Vector(self.ortho_matrix[:, 1].A.flatten(), p0=self.offset)
        return Vector(self.cell_matrix[1].A.flatten(), p0=self.offset)

    @property
    def a3(self):
        """Lattice vector :math:`\\mathbf{a}_3=\\mathbf{c}`."""
        # return Vector(self.ortho_matrix[:, 2].A.flatten(), p0=self.offset)
        return Vector(self.cell_matrix[2].A.flatten(), p0=self.offset)

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
                 cell_matrix=None, orientation_matrix=None, offset=None):

        direct_lattice = \
            Crystal3DLattice(a=a_star, b=b_star, c=c_star, alpha=alpha_star,
                             beta=beta_star, gamma=gamma_star,
                             a1=b1, a2=b2, a3=b3, cell_matrix=cell_matrix,
                             orientation_matrix=orientation_matrix)
        super().__init__(direct_lattice, nd=3, offset=offset)

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
