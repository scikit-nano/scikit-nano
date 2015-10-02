# -*- coding: utf-8 -*-
"""
=======================================================================
3D geometric regions (:mod:`sknano.core.geometric_regions._3D_regions`)
=======================================================================

.. currentmodule:: sknano.core.geometric_regions._3D_regions

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

import numbers
import numpy as np

from sknano.core.math import Point, Vector, vector as vec
from ._base import GeometricRegion, GeometricTransformsMixin

__all__ = ['Geometric3DRegion', 'Parallelepiped', 'Cuboid', 'Cube',
           'Ellipsoid', 'Sphere', 'Cylinder', 'Cone']


class Geometric3DRegion(GeometricRegion, GeometricTransformsMixin,
                        metaclass=ABCMeta):
    """Abstract base class for representing 3D geometric regions."""

    @property
    @abstractmethod
    def volume(self):
        """Volume of 3D geometric region."""
        raise NotImplementedError

    @property
    def measure(self):
        """Alias for :attr:`~Geometric3DRegion.volume`, which is the measure \
            of a 3D geometric region."""
        return self.volume


class Parallelepiped(Geometric3DRegion):
    """`Geometric3DRegion` for a parallelepiped.

    .. versionadded:: 0.3.0

    Represents a parallelepiped with origin :math:`(o_x, o_y, o_z)` and
    direction vectors :math:`\\mathbf{u}=(u_x, u_y, u_z)`,
    :math:`\\mathbf{v}=(v_x, v_y, v_z)`, and
    :math:`\\mathbf{w}=(w_x, w_y, w_z)`.

    Parameters
    ----------
    o : array_like, optional
        Parallelepiped origin. If `None`, it defaults to `o=[0, 0, 0]`.
    u, v, w : array_like, optional
        Parallelepiped direction vectors stemming from origin `o`.
        If `None`, then the default values are `u=[1, 0, 0]`,
        `v=[1, 1, 0]`, and `w=[0, 1, 1]`.

    Notes
    -----
    :class:`Parallelepiped` represents the bounded region
    :math:`\\left\\{o+\\lambda_1\\mathbf{u}+\\lambda_2\\mathbf{v}+
    \\lambda_3\\mathbf{w}\\in R^3|0\\le\\lambda_i\\le 1\\right\\}`,
    where :math:`\\mathbf{u}`, :math:`\\mathbf{v}`, and :math:`\\mathbf{w}`
    have to be linearly independent.

    Calling `Parallelepiped` with no parameters is equivalent to
    `Parallelepiped`\ `(o=[0, 0, 0], u=[1, 0, 0], v=[1, 1, 0], w=[0, 1, 1])`.


    """
    def __init__(self, o=None, u=None, v=None, w=None):

        self._o = Point()
        self._u = Vector()
        self._v = Vector()
        self._w = Vector()

        super().__init__()

        if o is None:
            o = [0, 0, 0]
        self.o = o

        if u is None:
            u = [1, 0, 0]
        self.u = u

        if v is None:
            v = [1, 1, 0]
        self.v = v

        if w is None:
            w = [0, 1, 1]
        self.w = w

        self.points.append(self.o)
        self.vectors.extend([self.u, self.v, self.w])

        self.fmtstr = "o={o!r}, u={u!r}, v={v!r}, w={w!r}"

    @property
    def o(self):
        """3D point coordinates :math:`(o_x, o_y, o_z)` of origin.

        Returns
        -------
        :class:`Point`
            3D :class:`Point` coordinates :math:`(o_x, o_y, o_z)` of origin.
        """
        return self._o

    @o.setter
    def o(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self.vectors.translate(Vector(p0=self.o, p=Point(value)))
        self._o[:] = Point(value)

    @property
    def u(self):
        """3D direction vector :math:`\\mathbf{u}=(u_x, u_y, u_z)`, with \
            origin :attr:`~Parallelepiped.o`"""
        return self._u

    @u.setter
    def u(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._u[:] = Vector(value, p0=self.o)

    @property
    def v(self):
        """3D direction vector :math:`\\mathbf{v}=(v_x, v_y, v_z)`, with \
            origin :attr:`~Parallelepiped.o`"""
        return self._v

    @v.setter
    def v(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._v[:] = Vector(value, p0=self.o)

    @property
    def w(self):
        """3D direction vector :math:`\\mathbf{w}=(w_x, w_y, w_z)`, with \
            origin :attr:`~Parallelepiped.o`"""
        return self._w

    @w.setter
    def w(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._w[:] = Vector(value, p0=self.o)

    @property
    def volume(self):
        """Parallelepiped volume, :math:`V`.

        Computed as:

        .. math::

           V = |\\mathbf{u}\\cdot\\mathbf{v}\\times\\mathbf{w}|

        """
        return np.abs(vec.scalar_triple_product(self.u, self.v, self.w))

    @property
    def centroid(self):
        """Parallelepiped centroid, :math:`(c_x, c_y, c_z)`.

        Computed as the 3D point :math:`(c_x, c_y, c_z)` with coordinates:

        .. math::

           c_x = o_x + \\frac{u_x + v_x + w_x}{2}

           c_y = o_y + \\frac{u_y + v_y + w_y}{2}

           c_z = o_z + \\frac{u_z + v_z + w_z}{2}

        where :math:`(o_x, o_y, o_z)`, :math:`(u_x, u_y, u_z)`,
        :math:`(v_x, v_y, v_z)`, and :math:`(w_x, w_y, w_z)`
        are the :math:`(x, y, z)` coordinates of the origin :math:`o`
        and :math:`(x, y, z)` components of the direction vectors
        :math:`\\mathbf{u}`, :math:`\\mathbf{v}`, and :math:`\\mathbf{w}`,
        respectively.

        Returns
        -------
        :class:`Point`
            3D :class:`Point` of centroid.

        """
        ox, oy, oz = self.o
        ux, uy, uz = self.u
        vx, vy, vz = self.v
        wx, wy, wz = self.w

        cx = ox + (ux + vx + wx) / 2
        cy = oy + (uy + vy + wy) / 2
        cz = oz + (uz + vz + wz) / 2

        return Point([cx, cy, cz])

    def contains(self, point):
        """Test region membership of `point` in :class:`Parallelepiped`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Parallelepiped`,
            `False` otherwise

        Notes
        -----
        A `point` :math:`(p_x, p_y, p_z)` is within the bounded region
        of a parallelepiped with origin :math:`(o_x, o_y, o_z)` and
        direction vectors :math:`\\mathbf{u}=(u_x, u_y, u_z)`,
        :math:`\\mathbf{v}=(v_x, v_y, v_z)`, and
        :math:`\\mathbf{w}=(w_x, w_y, w_z)` if the following is true:

        .. math::

           0\\le\\frac{
           v_z (w_x (p_y - o_y) + w_y (o_x - p_x)) +
           w_z (v_x (o_y - p_y) + v_y (p_x - o_x)) +
           o_z (v_y w_x - v_x w_y) + p_z (v_x w_y - v_y w_x)}{
           u_z (v_x w_y - v_y w_x) +
           u_y (v_z w_x - v_x w_z) +
           u_x (v_y w_z - v_z w_y)}\\le 1 \\land

           0\\le\\frac{
           u_z (w_x (p_y - o_y) + w_y (o_x - p_x)) +
           w_z (u_x (o_y - p_y) + u_y (p_x - o_x)) +
           o_z (u_y w_x - u_x w_y) +
           p_z (u_x w_y - u_y w_x)}{
           u_z (v_y w_x - v_x w_y) +
           u_y (v_x w_z - v_z w_x) +
           u_x (v_z w_y - v_y w_z)}\\le 1 \\land

           0\\le\\frac{
           u_z (v_x (p_y - o_y) + v_y (o_x - p_x)) +
           v_z (u_x (o_y - p_y) + u_y (p_x - o_x)) +
           o_z (u_y v_x - u_x v_y) +
           p_z (u_x v_y - u_y v_x)}{
           u_z (v_x w_y - v_y w_x) +
           u_y (v_z w_x - v_x w_z) +
           u_x (v_y w_z - v_z w_y)}\\le 1

        """
        px, py, pz = Point(point)

        ox, oy, oz = self.o
        ux, uy, uz = self.u
        vx, vy, vz = self.v
        wx, wy, wz = self.w

        q1 = (vz * (wx * (py - oy) + wy * (ox - px)) +
              wz * (vx * (oy - py) + vy * (px - ox)) +
              oz * (vy * wx - vx * wy) +
              pz * (vx * wy - vy * wx)) / \
            (uz * (vx * wy - vy * wx) +
             uy * (vz * wx - vx * wz) +
             ux * (vy * wz - vz * wy))

        q2 = (uz * (wx * (py - oy) + wy * (ox - px)) +
              wz * (ux * (oy - py) + uy * (px - ox)) +
              oz * (uy * wx - ux * wy) +
              pz * (ux * wy - uy * wx)) / \
            (uz * (vy * wx - vx * wy) +
             uy * (vx * wz - vz * wx) +
             ux * (vz * wy - vy * wz))

        q3 = (uz * (vx * (py - oy) + vy * (ox - px)) +
              vz * (ux * (oy - py) + uy * (px - ox)) +
              oz * (uy * vx - ux * vy) +
              pz * (ux * vy - uy * vx)) / \
            (uz * (vx * wy - vy * wx) +
             uy * (vz * wx - vx * wz) +
             ux * (vy * wz - vz * wy))

        return q1 >= 0 and q1 <= 1 and q2 >= 0 and q2 <= 1 and \
            q3 >= 0 and q3 <= 1

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Parallelepiped` \
            constructor parameters."""
        return dict(o=self.o, u=self.u, v=self.v, w=self.w)


class Cuboid(Geometric3DRegion):
    """`Geometric3DRegion` for a cuboid.

    .. versionadded:: 0.3.0

    Represents an axis-aligned cuboid with lower corner
    :math:`p_{\\mathrm{min}}=
    (x_{\\mathrm{min}},y_{\\mathrm{min}},z_{\\mathrm{min}})` and
    upper corner
    :math:`p_{\\mathrm{max}}=
    (x_{\\mathrm{max}},y_{\\mathrm{max}},z_{\\mathrm{max}})`.

    Parameters
    ----------
    pmin, pmax : array_like, optional
        The minimum and maximum 3D point coordinates of the axis-aligned
        cuboid from `pmin=[xmin, ymin, zmin]` to `pmax=[xmax, ymax, zmax]`.
    xmin, ymin, zmin : float, optional
        The minimum :math:`(x, y, z)` coordinates of the axis-aligned cuboid.
    xmax, ymax, zmax : float, optional
        The maximum :math:`(x, y, z)` coordinates of the axis-aligned cuboid.

    Notes
    -----
    :class:`Cuboid` represents the region
    :math:`\\left\\{\\{x, y, z\\}|
    x_{\\mathrm{min}}\\le x\\le x_{\\mathrm{max}}\\land
    y_{\\mathrm{min}}\\le y\\le y_{\\mathrm{max}}\\land
    z_{\\mathrm{min}}\\le z\\le z_{\\mathrm{max}}\\right\\}`

    Calling `Cuboid` with no parameters is equivalent to
    `Cuboid`\ `(pmin=[0, 0, 0], pmax=[1, 1, 1])`.

    """
    def __init__(self, pmin=None, pmax=None, xmin=0, ymin=0, zmin=0,
                 xmax=1, ymax=1, zmax=1):
        super().__init__()
        self._pmin = Point()
        self._pmax = Point()

        if pmin is None:
            pmin = [xmin, ymin, zmin]
        self.pmin = pmin

        if pmax is None:
            pmax = [xmax, ymax, zmax]
        self.pmax = pmax

        # assert self.pmin <= self.pmax
        assert np.all(np.less_equal(self.pmin, self.pmax))

        self.points.extend([self.pmin, self.pmax])
        self.fmtstr = "pmin={pmin!r}, pmax={pmax!r}"

    @property
    def pmin(self):
        """3D :class:`Point` at \
            (:attr:`~Cuboid.xmin`, :attr:`~Cuboid.ymin`, :attr:`~Cuboid.zmin`)
        """
        return self._pmin

    @pmin.setter
    def pmin(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._pmin[:] = Point(value)

    @property
    def pmax(self):
        """3D :class:`Point` at \
            (:attr:`~Cuboid.xmax`, :attr:`~Cuboid.ymax`, :attr:`~Cuboid.zmax`)
        """
        return self._pmax

    @pmax.setter
    def pmax(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._pmax[:] = Point(value)

    @property
    def xmin(self):
        """:math:`x_{\\mathrm{min}}` coordinate."""
        return self.pmin.x

    @xmin.setter
    def xmin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmin.x = float(value)

    @property
    def ymin(self):
        """:math:`y_{\\mathrm{min}}` coordinate."""
        return self.pmin.y

    @ymin.setter
    def ymin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmin.y = float(value)

    @property
    def zmin(self):
        """:math:`z_{\\mathrm{min}}` coordinate."""
        return self.pmin.z

    @zmin.setter
    def zmin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmin.z = float(value)

    @property
    def xmax(self):
        """:math:`x_{\\mathrm{max}}` coordinate."""
        return self.pmax.x

    @xmax.setter
    def xmax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmax.x = float(value)

    @property
    def ymax(self):
        """:math:`y_{\\mathrm{max}}` coordinate."""
        return self.pmax.y

    @ymax.setter
    def ymax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmax.y = float(value)

    @property
    def zmax(self):
        """:math:`z_{\\mathrm{max}}` coordinate."""
        return self.pmax.z

    @zmax.setter
    def zmax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmax.z = float(value)

    @property
    def a(self):
        """Distance between :math:`x_{\\mathrm{max}}-x_{\\mathrm{min}}`."""
        return self.xmax - self.xmin

    @property
    def b(self):
        """Distance between :math:`y_{\\mathrm{max}}-y_{\\mathrm{min}}`."""
        return self.ymax - self.ymin

    @property
    def c(self):
        """Distance between :math:`z_{\\mathrm{max}}-z_{\\mathrm{min}}`."""
        return self.zmax - self.zmin

    @property
    def volume(self):
        """:class:`Cuboid` volume, :math:`V=abc`."""
        return self.a * self.b * self.c

    @property
    def centroid(self):
        """:class:`Cuboid` centroid, :math:`(c_x, c_y, c_z)`.

        Computed as the 3D :class:`Point` :math:`(c_x, c_y, c_z)` with
        coordinates:

        .. math::

           c_x = \\frac{x_{\\mathrm{min}} + x_{\\mathrm{max}}}{2}

           c_y = \\frac{y_{\\mathrm{min}} + y_{\\mathrm{max}}}{2}

           c_z = \\frac{z_{\\mathrm{min}} + z_{\\mathrm{max}}}{2}

        Returns
        -------
        :class:`Point`
            3D :class:`Point` of centroid.

        """
        cx = (self.xmin + self.xmax) / 2
        cy = (self.ymin + self.ymax) / 2
        cz = (self.zmin + self.zmax) / 2
        return Point([cx, cy, cz])

    def contains(self, point):
        """Test region membership of `point` in :class:`Cuboid`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Cuboid`,
            `False`, otherwise.

        Notes
        -----
        A point :math:`(p_x, p_y, p_z)` is within the bounded region of a
        cuboid with lower corner at
        :math:`p_{\\mathrm{min}}=
        (x_{\\mathrm{min}}, y_{\\mathrm{min}}, z_{\\mathrm{min}})` and
        upper corner at
        :math:`p_{\\mathrm{max}}=
        (x_{\\mathrm{max}}, y_{\\mathrm{max}}, z_{\\mathrm{max}})` if the
        following is true:

        .. math::

           x_{\mathrm{min}}\\le x\\le x_{\\mathrm{max}}\\land

           y_{\mathrm{min}}\\le y\\le y_{\\mathrm{max}}\\land

           z_{\mathrm{min}}\\le z\\le z_{\\mathrm{max}}

        """
        px, py, pz = Point(point)

        return (px >= self.xmin) and (px <= self.xmax) and \
            (py >= self.ymin) and (py <= self.ymax) and \
            (pz >= self.zmin) and (pz <= self.zmax)

    def rotate(self, **kwargs):
        super().rotate(**kwargs)
        pmin = self.pmin.copy()
        pmax = self.pmax.copy()
        self.pmin = np.minimum(pmin, pmax)
        self.pmax = np.maximum(pmin, pmax)
        assert np.all(np.less_equal(self.pmin, self.pmax))

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Cuboid` \
            constructor parameters."""
        return dict(pmin=self.pmin, pmax=self.pmax)


class Cube(Geometric3DRegion):
    """`Geometric3DRegion` for a cube.

    .. versionadded:: 0.3.0

    :class:`Cube` represents the region
    :math:`\\left\\{\\{c_i\\pm\\frac{a}{2}\\}|a>0\\forall
    i\\in\\{x,y,z\\}\\right\\}`.

    Parameters
    ----------
    center : array_like, optional
        The :math:`(x,y,z)` coordinate of the axis-aligned cube
        center :math:`(c_x, c_y, c_z)`.
    a : float, optional
        Side length :math:`a` of axis-aligned cube.

    Notes
    -----
    Calling :class:`Cube` with no parameters is equivalent to
    :class:`Cube`\ `(center=[0, 0, 0], a=1)`.

    """
    def __init__(self, center=None, a=1):
        super().__init__()

        self._center = Point()

        if center is None:
            center = [0, 0, 0]
        self.center = center
        self.a = a
        self.points.append(self.center)

        self.fmtstr = "center={center!r}, a={a:.3f}"

    @property
    def center(self):
        """Center point :math:`(c_x, c_y, c_z)` of axis-aligned cube."""
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._center[:] = Point(value)

    @property
    def a(self):
        """Side length :math:`a` of the axis-aligned cube."""
        return self._a

    @a.setter
    def a(self, value):
        self._a = float(value)

    @property
    def centroid(self):
        """Alias for :attr:`~Cube.center`."""
        return self.center

    @property
    def volume(self):
        """:class:`Cube` volume, :math:`V=a^3`."""
        return self.a ** 3

    def contains(self, point):
        """Test region membership of `point` in :class:`Cube`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Cube`,
            `False` otherwise

        Notes
        -----
        A `point` :math:`(p_x, p_y, p_z)` is within the bounded region of
        a cube with center :math:`(c_x, c_y, c_z)` and side length
        :math:`a` if the following is true:

        .. math::

           c_i-\\frac{a}{2}\\le p_i\\le c_i+\\frac{a}{2}\\forall
           i\\in \\{x, y, z\\}

        """
        px, py, pz = Point(point)
        cx, cy, cz = self.center
        a = self.a
        xmin = cx - a / 2
        xmax = cx + a / 2
        ymin = cy - a / 2
        ymax = cy + a / 2
        zmin = cz - a / 2
        zmax = cz + a / 2
        return (px >= xmin) and (px <= xmax) and \
            (py >= ymin) and (py <= ymax) and \
            (pz >= zmin) and (pz <= zmax)

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Cube` \
            constructor parameters."""
        return dict(center=self.center, a=self.a)


class Ellipsoid(Geometric3DRegion):
    """`Geometric3DRegion` for an ellipsoid.

    .. versionadded:: 0.3.0

    Represents an axis-aligned ellipsoid centered at the point
    :math:`(c_x, c_y, c_z)` with semi-axes lengths
    :math:`r_x, r_y, r_z`.

    Parameters
    ----------
    center : array_like or :class:`~sknano.core.Point`
        3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the center point
        coordinates :math:`(x,y,z)` of the axis-aligned
        :class:`Ellipsoid` with semi-axes lengths `rx`, `ry`, `rz`.
    rx, ry, rz : float, optional
        Semi-axes lengths :math:`r_x, r_y, r_z` of axis-aligned
        :class:`Ellipsoid` centered at `center`.

    Notes
    -----
    :class:`Ellipsoid` represents the axis-aligned ellipsoid:

    .. math::

       \\left\\{\\{x, y, z\\}\\in R^3|
       \\left(\\frac{x-c_x}{r_x}\\right)^2 +
       \\left(\\frac{y-c_y}{r_y}\\right)^2 +
       \\left(\\frac{z-c_z}{r_z}\\right)^2\\le 1\\right\\}

    Calling :class:`Ellipsoid` with no parameters is equivalent to
    :class:`Ellipsoid`\ `(center=[0, 0, 0], rx=1, ry=1, rz=1)`.

    """
    def __init__(self, center=None, rx=1, ry=1, rz=1):
        super().__init__()

        self._center = Point()

        if center is None:
            center = [0, 0, 0]
        self.center = center

        self.rx = rx
        self.ry = ry
        self.rz = rz

        self.points.append(self.center)

        self.fmtstr = \
            "center={center!r}, rx={rx:.3f}, ry={ry:.3f}, rz={rz:.3f}"

    @property
    def center(self):
        """:class:`Ellipsoid` center point :math:`(c_x, c_y, c_z)`."""
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._center[:] = Point(value)

    @property
    def rx(self):
        """Length of semi-axis :math:`r_x`."""
        return self._rx

    @rx.setter
    def rx(self, value):
        self._rx = float(value)

    @property
    def ry(self):
        """Length of semi-axis :math:`r_y`."""
        return self._ry

    @ry.setter
    def ry(self, value):
        self._ry = float(value)

    @property
    def rz(self):
        """Length of semi-axis :math:`r_z`."""
        return self._rz

    @rz.setter
    def rz(self, value):
        self._rz = float(value)

    @property
    def centroid(self):
        """Alias for :attr:`~Ellipsoid.center`."""
        return self.center

    @property
    def volume(self):
        """:class:`Ellipsoid` volume, \
            :math:`V=\\frac{4}{3}\\pi r_x r_y r_z`."""
        return 4 / 3 * np.pi * self.rx * self.ry * self.rz

    def contains(self, point):
        """Test region membership of `point` in :class:`Ellipsoid`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Ellipsoid`,
            `False` otherwise

        Notes
        -----
        A `point` :math:`(p_x, p_y, p_z)` is within the bounded region of
        an ellipsoid with center :math:`(c_x, c_y, c_z)` and semi-axes lengths
        :math:`r_x, r_y, r_z` if the following is true:

        .. math::

           \\left(\\frac{p_x - c_x}{r_x}\\right)^2 +
           \\left(\\frac{p_y - c_y}{r_y}\\right)^2 +
           \\left(\\frac{p_z - c_z}{r_z}\\right)^2\\le 1

        """
        px, py, pz = Point(point)
        cx, cy, cz = self.center
        rx, ry, rz = self.rx, self.ry, self.rz

        q1 = (px - cx) ** 2 / rx ** 2
        q2 = (py - cy) ** 2 / ry ** 2
        q3 = (pz - cz) ** 2 / rz ** 2

        return q1 + q2 + q3 <= 1.0

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Ellipsoid` \
            constructor parameters."""
        return dict(center=self.center, rx=self.rx, ry=self.ry, rz=self.rz)


class Sphere(Geometric3DRegion):
    """`Geometric3DRegion` for a sphere.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : array_like, optional
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Sphere` center.
    r : float, optional
        Sphere radius :math:`r`
    """
    def __init__(self, center=None, r=1):
        super().__init__()

        self._center = Point()

        if center is None:
            center = [0, 0, 0]
        self.center = center

        self.r = r

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, r={r:.3f}"

    @property
    def center(self):
        """:class:`Sphere` center point :math:`(h, k, l)`."""
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._center[:] = Point(value)

    @property
    def centroid(self):
        """Alias for :attr:`~Sphere.center`."""
        return self.center

    @property
    def r(self):
        """:class:`Sphere` radius, :math:`r`."""
        return self._r

    @r.setter
    def r(self, value):
        self._r = float(value)

    @property
    def volume(self):
        """:class:`Sphere` volume, :math:`V=\\frac{4}{3}\\pi r^3`."""
        return 4 / 3 * np.pi * self.r ** 3

    def contains(self, point):
        """Test region membership of `point` in :class:`Sphere`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Sphere`,
            `False` otherwise.

        Notes
        -----
        A `point` :math:`(p_x, p_y, p_z)` is within the bounded region
        of a sphere with center :math:`(h, k, l)` and radius :math:`r`
        if the following is true:

        .. math::

           (p_x - h)^2 + (p_y - k)^2 + (p_z - l)^2 \\le r^2

        """
        x, y, z = Point(point)
        h, k, l = self.center
        r = self.r

        return (x - h) ** 2 + (y - k) ** 2 + (z - l) ** 2 <= r ** 2

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Sphere` \
            constructor parameters."""
        return dict(center=self.center, r=self.r)


class Cylinder(Geometric3DRegion):
    """`Geometric3DRegion` for a cylinder.

    .. versionadded:: 0.3.10

    Represents a cylinder of radius :math:`r` around the line from
    :math:`(x_1, y_1, z_1)` to :math:`(x_2, y_2, z_2)`.

    Parameters
    ----------
    p1, p2 : array_like, optional
        3-tuples or :class:`~sknano.core.Point` class instances
        specifying the `Cylinder` axis from point :math:`p_1=(x_1,y_1,z_1)`
        to :math:`p_2=(x_2,y_2,z_2)`.
    r : float, optional
        `Cylinder` radius :math:`r`

    Notes
    -----
    :class:`Cylinder` represents a cylinder region
    :math:`\\left\\{p_1+\\rho\\cos(\\theta)\\mathbf{v}_1 +
    \\rho\\sin(\\theta)\\mathbf{v}_2 + \\mathbf{v}_3 z|
    0\\le\\theta\\le 2\\pi\\land
    0\\le\\rho\\le 1\\land
    0\\le z\\le 1\\right\\}`
    where :math:`\\mathbf{v}_3=p_2 - p_1` and the vectors
    :math:`\\{\\mathbf{v}_1,\\mathbf{v}_2,\\mathbf{v}_3\\}` are orthogonal
    with :math:`|\\mathbf{v}_1|=|\\mathbf{v}_2|=1`, and
    :math:`p_1=(x_1, y_1, z_1)` and :math:`p_2=(x_2,y_2,z_2)`.

    Calling :class:`Cylinder` with no parameters is equivalent to
    :class:`Cylinder`\ `(p1=[0, 0, -1], p2=[0, 0, 1], r=1)`.

    """
    def __init__(self, p1=None, p2=None, r=1):
        super().__init__()

        self._p1 = Point()
        self._p2 = Point()

        if p1 is None:
            p1 = [0, 0, -1]
        if p2 is None:
            p2 = [0, 0, 1]

        self.p1 = p1
        self.p2 = p2
        self.points.extend([self.p1, self.p2])
        self.r = r
        self.fmtstr = "p1={p1!r}, p2={p2!r}, r={r:.3f}"

    @property
    def r(self):
        """:class:`Cylinder` radius :math:`r`."""
        return self._r

    @r.setter
    def r(self, value):
        self._r = float(value)

    @property
    def p1(self):
        """:class:`Cylinder` axis point :math:`p_1=(x_1, y_1, z_1)`."""
        return self._p1

    @p1.setter
    def p1(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._p1[:] = Point(value)

    @property
    def p2(self):
        """:class:`Cylinder` axis point :math:`p_2=(x_2, y_2, z_2)`."""
        return self._p2

    @p2.setter
    def p2(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._p2[:] = Point(value)

    @property
    def axis(self):
        """:class:`Cylinder` axis :class:`Vector` \
            :math:`\\boldsymbol{\\ell}=p_2 - p_1`.

        Returns
        -------
        :class:`Vector`
            3D :class:`Vector` along :class:`Cylinder` axis from
            :class:`Point` :math:`p_1=(x_1, y_1, z_1)` to
            :class:`Point` :math:`p_2=(x_2, y_2, z_2)`.

        """
        return Vector(p0=self.p1, p=self.p2)

    @property
    def centroid(self):
        """:class:`Cylinder` centroid, :math:`(c_x, c_y, c_z)`.

        Computed as:

        .. math::

           c_x = \\frac{x_1 + x_2}{2}

           c_y = \\frac{y_1 + y_2}{2}

           c_z = \\frac{z_1 + z_2}{2}

        """
        return (self.p1 + self.p2) / 2

    @property
    def volume(self):
        """:class:`Cylinder` volume, :math:`V=\\pi r^2 \\ell`."""
        return np.pi * self.r ** 2 * self.axis.length

    def contains(self, point):
        """Test region membership of `point` in :class:`Cylinder`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Cylinder`,
            `False` otherwise.

        Notes
        -----
        A `point` :math:`(p_x, p_y, p_z)` is within the bounded region
        of a cylinder with radius :math:`r` around the line from
        :math:`p_1=(x_1, y_1, z_1)` to :math:`p_2 = (x_2, y_2, z_2)`
        if the following is true:

        .. math::

           0\\le q\\le 1\\land
           (x_1 - p_x + (x_2 - x_1) q)^2 +
           (y_1 - p_y + (y_2 - y_1) q)^2 +
           (z_1 - p_z + (z_2 - z_1) q)^2 \\le r^2

        where :math:`q` is:

        .. math::

           q = \\frac{(p_x - x_1)(x_2 - x_1) +
           (p_y - y_1)(y_2 - y_1) +
           (p_z - z_1)(z_2 - z_1)}{(x_2 - x_1)^2 +
           (y_2 - y_1)^2 + (z_2 - z_1)^2}

        """
        px, py, pz = Point(point)
        x1, y1, z1 = self.p1
        x2, y2, z2 = self.p2
        r = self.r

        q1 = ((px - x1) * (x2 - x1) +
              (py - y1) * (y2 - y1) +
              (pz - z1) * (z2 - z1)) / \
            ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

        q2 = (x1 - px + (x2 - x1) * q1) ** 2 + \
            (y1 - py + (y2 - y1) * q1) ** 2 + \
            (z1 - pz + (z2 - z1) * q1) ** 2

        return (not np.allclose(x1, x2) or not np.allclose(y1, y2) or
                not np.allclose(z1, z2)) and r > 0 and q1 >= 0 and q1 <= 1 \
            and q2 <= r ** 2

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Cylinder` \
            constructor parameters."""
        return dict(p1=self.p1, p2=self.p2, r=self.r)


class Cone(Geometric3DRegion):
    """`Geometric3DRegion` for a cone.

    .. versionadded:: 0.3.10

    Represents a cone with a base of radius :math:`r` centered at
    :math:`p_1=(x_1,y_1,z_1)` and a tip at :math:`p_2=(x_2, y_2, z_2)`.

    Parameters
    ----------
    p1, p2 : array_like, optional
        3-tuples or :class:`~sknano.core.Point` class instances
        for a :class:`Cone` with a base of radius `r` centered at
        :math:`p_1=(x_1,y_1,z_1)` and a tip at :math:`p_2=(x_2,y_2,z_2)`.
    r : float, optional
        Radius :math:`r` of :class:`Cone` base

    Notes
    -----
    :class:`Cone` represents a bounded cone region
    :math:`\\left\\{p_1+\\rho(1-z)\\cos(\\theta)\\mathbf{v}_1 +
    \\rho(1-z)\\sin(\\theta)\\mathbf{v}_2+\\mathbf{v}_3 z|
    0\\le\\theta\\le 2\pi\\land
    0\\le\\rho\\le 1\\land
    0\\le z\\le 1\\right\\}`
    where :math:`\\mathbf{v}_3=p_2-p_1` and vectors
    :math:`(\\mathbf{v}_1,\\mathbf{v}_2,\\mathbf{v}_3)` are orthogonal
    with :math:`|\\mathbf{v}_1|=|\\mathbf{v}_2|=1` and
    :math:`p_1=(x_1,y_1,z_1)` and :math:`p_2=(x_2, y_2, z_2)`.

    Calling :class:`Cone` with no parameters is equivalent to
    :class:`Cone`\ `(p1=[0, 0, 0], p2=[0, 0, 2], r=1)`.

    """
    def __init__(self, p1=None, p2=None, r=1):
        super().__init__()

        self._p1 = Point()
        self._p2 = Point()

        if p1 is None:
            p1 = [0, 0, 0]
        if p2 is None:
            p2 = [0, 0, 2]

        self.p1 = p1
        self.p2 = p2
        self.points.extend([self.p1, self.p2])
        self.r = r
        self.fmtstr = "p1={p1!r}, p2={p2!r}, r={r:.3f}"

    @property
    def r(self):
        """Radius :math:`r` of :class:`Cone` base."""
        return self._r

    @r.setter
    def r(self, value):
        self._r = float(value)

    @property
    def p1(self):
        """Center point :math:`(x_1, y_1, z_1)` of :class:`Cone` base."""
        return self._p1

    @p1.setter
    def p1(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._p1[:] = Point(value)

    @property
    def p2(self):
        """Point :math:`(x_2, y_2, z_2)` of :class:`Cone` tip."""
        return self._p2

    @p2.setter
    def p2(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._p2[:] = Point(value)

    @property
    def axis(self):
        """:class:`Cone` axis :class:`Vector` from \
            :math:`\\boldsymbol{\\ell}=p_2 - p_1`

        Returns
        -------
        :class:`Vector`

        """
        return Vector(p0=self.p1, p=self.p2)

    @property
    def centroid(self):
        """:class:`Cone` centroid, :math:`(c_x, c_y, c_z)`.

        Computed as:

        .. math::

           c_x = \\frac{3 x_1 + x_2}{4}

           c_y = \\frac{3 y_1 + y_2}{4}

           c_z = \\frac{3 z_1 + z_2}{4}

        Returns
        -------
        :class:`Point`
            3D :class:`Point` of :class:`Cone` centroid.

        """
        return (3 * self.p1 + self.p2) / 4

    @property
    def volume(self):
        """:class:`Cone` volume, :math:`V=\\frac{1}{3}\\pi r^2 \\ell`."""
        return np.pi * self.r ** 2 * self.axis.length / 3

    def contains(self, point):
        """Test region membership of `point` in :class:`Cone`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Cone`,
            `False` otherwise.

        Notes
        -----
        A `point` :math:`(p_x, p_y, p_z)` is within the bounded region
        of a cone with a base of radius :math:`r` centered at
        :math:`p_1=(x_1, y_1, z_1)` and tip at :math:`p_2 = (x_2, y_2, z_2)`
        if the following is true:

        .. math::

           0\\le q\\le 1\\land
           (x_1 - p_x + (x_2 - x_1) q)^2 +
           (y_1 - p_y + (y_2 - y_1) q)^2 +
           (z_1 - p_z + (z_2 - z_1) q)^2 \\le r^2 q^2

        where :math:`q` is:

        .. math::

           q = \\frac{(p_x - x_1)(x_2 - x_1) +
           (p_y - y_1)(y_2 - y_1) +
           (p_z - z_1)(z_2 - z_1)}{(x_2 - x_1)^2 +
           (y_2 - y_1)^2 + (z_2 - z_1)^2}

        """
        px, py, pz = Point(point)
        x1, y1, z1 = self.p1
        x2, y2, z2 = self.p2
        r = self.r

        q1 = ((px - x1) * (x2 - x1) +
              (py - y1) * (y2 - y1) +
              (pz - z1) * (z2 - z1)) / \
            ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

        q2 = (x1 - px + (x2 - x1) * q1) ** 2 + \
            (y1 - py + (y2 - y1) * q1) ** 2 + \
            (z1 - pz + (z2 - z1) * q1) ** 2

        q3 = r ** 2 * q1 ** 2

        return (not np.allclose(x1, x2) or not np.allclose(y1, y2) or
                not np.allclose(z1, z2)) and r > 0 and q1 >= 0 and q1 <= 1 \
            and q2 <= q3

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Cone` \
            constructor parameters."""
        return dict(p1=self.p1, p2=self.p2, r=self.r)
