# -*- coding: utf-8 -*-
"""
=======================================================================
3D geometric shapes (:mod:`sknano.utils.geometric_shapes._3D_shapes`)
=======================================================================

.. currentmodule:: sknano.utils.geometric_shapes._3D_shapes

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
           'Ellipsoid', 'Spheroid', 'Sphere']


class Geometric3DRegion(GeometricRegion, GeometricTransformsMixin,
                        metaclass=ABCMeta):
    """Abstract base class for representing 3D geometric regions."""

    @property
    @abstractmethod
    def volume(self):
        """Volume of 3D geometric region."""
        raise NotImplementedError


class Parallelepiped(Geometric3DRegion):
    """Abstract representation of parallelepiped.

    .. versionadded:: 0.3.0

    Represents a parallelepiped with origin :math:`o` and directions
    :math:`u, v, w`.

    Parameters
    ----------
    o : array_like, optional
        parallelepiped origin
    u, v, w : array_like, optional
        parallelepiped direction vectors stemming from origin `o`.


    """
    def __init__(self, o=None, u=None, v=None, w=None):
        super().__init__()

        if o is None:
            o = Point(nd=3)
        elif isinstance(o, (tuple, list, np.ndarray)):
            o = Point(o)
        self._o = o

        if u is None:
            u = [1., 0., 0.]
        self._u = Vector(u, p0=self._o)

        if v is None:
            v = [1., 1., 0.]
        self._v = Vector(v, p0=self._o)

        if w is None:
            w = [0., 1., 1.]
        self._w = Vector(w, p0=self._o)

        self.points.append(self.o)
        self.vectors.extend([self.u, self.v, self.w])

        self.fmtstr = "o={o!r}, u={u!r}, v={v!r}, w={w!r}"

    @property
    def o(self):
        return self._o

    @o.setter
    def o(self, value):
        self._o[:] = value

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        self._u[:] = value

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        self._v[:] = value

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, value):
        self._w[:] = value

    @property
    def center(self):
        return self.centroid

    @property
    def volume(self):
        u = self.u
        v = self.v
        w = self.w
        return vec.scalar_triple_product(u, v, w)

    @property
    def centroid(self):
        ox, oy, oz = self.o
        ux, uy, uz = self.u
        vx, vy, vz = self.v
        wx, wy, wz = self.w

        xcom = 0.5 * (2 * ox + ux + vx + wx)
        ycom = 0.5 * (2 * oy + uy + vy + wy)
        zcom = 0.5 * (2 * oy + uy + vy + wz)

        return Point([xcom, ycom, zcom])

    def contains(self, point):
        """Check if point is contained within volume of `Parallelepiped`."""
        x, y, z = Point(point)

        ox, oy, oz = self.o
        ux, uy, uz = self.u
        vx, vy, vz = self.v
        wx, wy, wz = self.w

        d1 = (vz * (wx * (y - oy) + wy * (ox - x)) +
              wz * (vx * (oy - y) + vy * (x - ox)) +
              oz * (vy * wx - vx * wy) +
              z * (vx * wy - vy * wx)) / \
            (uz * (vx * wy - vy * wx) +
             uy * (vz * wx - vx * wz) +
             ux * (vy * wz - vz * wy))

        d2 = (uz * (wx * (y - oy) + wy * (ox - x)) +
              wz * (ux * (oy - y) + uy * (x - ox)) +
              oz * (uy * wx - ux * wy) +
              z * (ux * wy - uy * wx)) / \
            (uz * (vy * wx - vx * wy) +
             uy * (vx * wz - vz * wx) +
             ux * (vz * wy - vy * wz))

        d3 = (uz * (vx * (y - oy) + vy * (ox - x)) +
              vz * (ux * (oy - y) + uy * (x - ox)) +
              oz * (uy * vx - ux * vy) +
              z * (ux * vy - uy * vx)) / \
            (uz * (vx * wy - vy * wx) +
             uy * (vz * wx - vx * wz) +
             ux * (vy * wz - vz * wy))

        return d1 >= 0 and d1 <= 1 and d2 >= 0 and d2 <= 1 and \
            d3 >= 0 and d3 <= 1

    def todict(self):
        return dict(o=self.o, u=self.u, v=self.v, w=self.w)


class Cuboid(Geometric3DRegion):
    """Abstract representation of cuboid.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    pmin, pmax : sequence, optional
    xmin, ymin, zmin : float, optional
    xmax, ymax, zmax : float, optional

    """
    def __init__(self, pmin=None, pmax=None,
                 xmin=None, ymin=None, zmin=None,
                 xmax=None, ymax=None, zmax=None):
        super().__init__()

        if pmin is None:
            pmin = Point([xmin, ymin, zmin])
        elif isinstance(pmin, (tuple, list, np.ndarray)):
            pmin = Point(pmin)
        self._pmin = pmin

        if pmax is None:
            pmax = Point([xmax, ymax, zmax])
        elif isinstance(pmax, (tuple, list, np.ndarray)):
            pmax = Point(pmax)
        self._pmax = pmax

        self.points.extend([self.pmin, self.pmax])
        self.fmtstr = "pmin={pmin!r}, pmax={pmax!r}"

    @property
    def pmin(self):
        return self._pmin

    @pmin.setter
    def pmin(self, point):
        if not isinstance(point, (list, np.ndarray)) or \
                len(point) != self.pmin.nd:
            raise TypeError('Expected a 3-element list or ndarray')
        self._pmin[:] = point

    @property
    def pmax(self):
        return self._pmax

    @pmax.setter
    def pmax(self, point):
        if not isinstance(point, (list, np.ndarray)) or \
                len(point) != self.pmax.nd:
            raise TypeError('Expected a 3-element list or ndarray')
        self._pmax[:] = point

    @property
    def xmin(self):
        return self.pmin.x

    @xmin.setter
    def xmin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pmin.x = value

    @property
    def ymin(self):
        return self.pmin.y

    @ymin.setter
    def ymin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pmin.y = value

    @property
    def zmin(self):
        return self.pmin.z

    @zmin.setter
    def zmin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pmin.z = value

    @property
    def xmax(self):
        return self.pmax.x

    @xmax.setter
    def xmax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pmax.x = value

    @property
    def ymax(self):
        return self.pmax.y

    @ymax.setter
    def ymax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pmax.y = value

    @property
    def zmax(self):
        return self.pmax.z

    @zmax.setter
    def zmax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pmax.z = value

    @property
    def center(self):
        return self.centroid

    @property
    def centroid(self):
        xcom = 0.5 * (self.xmin + self.xmax)
        ycom = 0.5 * (self.ymin + self.ymax)
        zcom = 0.5 * (self.zmin + self.zmax)
        return Point([xcom, ycom, zcom])

    @property
    def a(self):
        return self.xmax - self.xmin

    @property
    def b(self):
        return self.ymax - self.ymin

    @property
    def c(self):
        return self.zmax - self.zmin

    @property
    def volume(self):
        a = self.a
        b = self.b
        c = self.c
        return a * b * c

    def contains(self, point):
        """Check if point is contained within volume of `Cuboid`."""
        x, y, z = Point(point)

        return (x >= self.xmin) and (x <= self.xmax) and \
            (y >= self.ymin) and (y <= self.ymax) and \
            (z >= self.zmin) and (z <= self.zmax)

    def todict(self):
        return dict(pmin=self.pmin, pmax=self.pmax)


class Cube(Geometric3DRegion):
    """Abstract data structure representing a cube.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the cube center.
    a : float, optional
        length of edge

    """
    def __init__(self, center=None, a=None):
        super().__init__()

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if a is None:
            a = 1.0
        self._a = a

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, a={a:.3f}"

    @property
    def center(self):
        return self._center

    @property
    def a(self):
        return self._a

    @property
    def centroid(self):
        self.center

    @property
    def volume(self):
        a = self.a
        return a**3

    def contains(self, point):
        x, y, z = Point(point)
        h, k, l = self.center
        a = self.a
        xmin = h - a / 2
        ymin = k - a / 2
        zmin = l - a / 2
        xmax = h + a / 2
        ymax = k + a / 2
        zmax = l + a / 2
        return (x >= xmin) and (x <= xmax) and \
            (y >= ymin) and (y <= ymax) and \
            (z >= zmin) and (z <= zmax)

    def todict(self):
        return dict(center=self.center, a=self.a)


class Ellipsoid(Geometric3DRegion):
    """Abstract data structure representing an ellipsoid.

    .. versionadded:: 0.3.0

    The general ellipsoid is a quadratic surface with is given in
    Cartesian coordinates by:

    .. math::

       \\frac{x^2}{a^2} + \\frac{y^2}{b^2} + \\frac{z^2}{c^2} = 1

    Parameters
    ----------
    center : sequence or :class:`~sknano.core.Point`
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Ellipsoid` center.
    a, b, c : float, optional
        Semi-principal axes :math:`a, b, c` of axis-aligned
        :class:`Ellipsoid`

    """
    def __init__(self, center=None, a=None, b=None, c=None):
        super().__init__()

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if a is None:
            a = 0.5
        if b is None:
            b = 0.5
        if c is None:
            c = 0.5
        self._a, self._b, self._c = a, b, c

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, a={a:.3f}, b={b:.3f}, c={c:.3f}"

    @property
    def center(self):
        return self._center

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def centroid(self):
        return self.center

    @property
    def volume(self):
        a = self.a
        b = self.b
        c = self.c
        return 4 / 3 * np.pi * a * b * c

    def contains(self, point):
        """Check if point is contained within volume of :class:`Ellipsoid`."""
        x, y, z = Point(point)
        h, k, l = self.center
        a, b, c = self.a, self.b, self.c

        return (x - h)**2 / a**2 + (y - k)**2 / b**2 + \
            (z - l)**2 / c**2 <= 1.0

    def todict(self):
        return dict(center=self.center, a=self.a, b=self.b, c=self.c)


class Spheroid(Geometric3DRegion):
    """Abstract data structure representing a spheroid.

    .. versionadded:: 0.3.0

    The general spheroid is a quadratic surface with is given in
    Cartesian coordinates by:

    .. math::

       \\frac{x^2 + y^2}{a^2} + \\frac{z^2}{c^2} = 1

    Parameters
    ----------
    center : sequence, optional
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Spheroid` center.
    a, c : float, optional
        Semi-axes :math:`a, c` of axis-aligned
        :class:`Spheroid` with symmetry axis along :math:`z` axis.
    """
    def __init__(self, center=None, a=None, c=None):
        super().__init__()

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if a is None:
            a = 0.5
        if c is None:
            c = 0.5

        self._a, self._c = a, c

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, a={a:.3f}, c={c:.3f}"

    @property
    def center(self):
        return self._center

    @property
    def a(self):
        return self._a

    @property
    def c(self):
        return self._c

    @property
    def centroid(self):
        return self.center

    @property
    def volume(self):
        a = self.a
        c = self.c
        return 4 / 3 * np.pi * a**2 * c

    def contains(self, point=None):
        """Check if point is contained within volume of :class:`Spheroid`."""
        x, y, z = point

        h, k, l = self.center
        a, c = self.a, self.c

        return ((x - h)**2 + (y - k)**2) / a**2 + (z - l)**2 / c**2 <= 1.0

    def todict(self):
        return dict(center=self.center, a=self.a, c=self.c)


class Sphere(Geometric3DRegion):
    """Abstract data structure representing a sphere.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Sphere` center.
    r : float, optional
        Sphere radius :math:`r`
    """
    def __init__(self, center=None, r=None):
        super().__init__()

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if r is None:
            r = 1.0
        self._r = r

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, r={r:.3f}"

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value

    @property
    def center(self):
        return self._center

    @property
    def centroid(self):
        return self.center

    @property
    def volume(self):
        r = self.r
        return 4 / 3 * np.pi * r**3

    def contains(self, point):
        x, y, z = point
        h, k, l = self.center
        r = self.r

        return (x - h)**2 + (y - k)**2 + (z - l)**2 <= r**2

    def todict(self):
        return dict(center=self.center, r=self.r)
