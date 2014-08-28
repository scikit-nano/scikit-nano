# -*- coding: utf-8 -*-
"""
=======================================================================
3D geometric shapes (:mod:`sknano.utils.geometric_shapes._3D_shapes`)
=======================================================================

.. currentmodule:: sknano.utils.geometric_shapes._3D_shapes

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Point, Vector, vector as vec
from ._base import GeometricRegion

__all__ = ['Geometric3DRegion', 'Parallelepiped', 'Cuboid', 'Cube',
           'Ellipsoid', 'Spheroid', 'Sphere']


class Geometric3DRegion(GeometricRegion):
    """Abstract base class for representing 3D geometric regions."""
    __metaclass__ = ABCMeta

    @abstractproperty
    def volume(self):
        """Volume of 3D geometric region."""
        raise NotImplementedError


class GeometricTransformsMixin(
    def rotate(self):
        pass

    def translate(self):
        pass


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

        if o is None:
            o = Point(nd=3)
        elif isinstance(o, (tuple, list, np.ndarray)):
            o = Point(o)
        self._o = o

        if u is None:
            u = Vector([1., 0., 0.])
        elif isinstance(u, (tuple, list, np.ndarray)):
            u = Vector(u)
        self._u = u

        if v is None:
            v = Vector([0., 1., 0.])
        elif isinstance(v, (tuple, list, np.ndarray)):
            v = Vector(v)
        self._v = v

        if w is None:
            w = Vector([1., 1., 1.])
        elif isinstance(w, (tuple, list, np.ndarray)):
            w = Vector(w)
        self._w = w

    def __repr__(self):
        return "Parallelepiped(o={!r}, u={!r}, v={!r}, w={!r})".format(
            self.o.tolist(), self.u.tolist(), self.v.tolist(), self.w.tolist())

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
        o = self.o
        u = self.u
        v = self.v
        w = self.w

        xcom = 0.5 * (2 * o.x + u.x + v.x + w.x)
        ycom = 0.5 * (2 * o.y + u.y + v.y + w.y)
        zcom = 0.5 * (2 * o.y + u.y + v.y + w.z)

        return Point([xcom, ycom, zcom])

    def contains_point(self, point):
        """Check if point is contained within volume of `Parallelepiped`."""
        p = Point(point)

        o = self.o
        u = self.u
        v = self.v
        w = self.w

        d1 = (v.z * (w.x * (p.y - o.y) + w.y * (o.x - p.x)) +
              w.z * (v.x * (o.y - p.y) + v.y * (p.x - o.x)) +
              o.z * (v.y * w.x - v.x * w.y) +
              p.z * (v.x * w.y - v.y * w.x)) / \
            (u.z * (v.x * w.y - v.y * w.x) +
             u.y * (v.z * w.x - v.x * w.z) +
             u.x * (v.y * w.z - v.z * w.y))

        d2 = (u.z * (w.x * (p.y - o.y) + w.y * (o.x - p.x)) +
              w.z * (u.x * (o.y - p.y) + u.y * (p.x - o.x)) +
              o.z * (u.y * w.x - u.x * w.y) +
              p.z * (u.x * w.y - u.y * w.x)) / \
            (u.z * (v.y * w.x - v.x * w.y) +
             u.y * (v.x * w.z - v.z * w.x) +
             u.x * (v.z * w.y - v.y * w.z))

        d3 = (u.z * (v.x * (p.y - o.y) + v.y * (o.x - p.x)) +
              v.z * (u.x * (o.y - p.y) + u.y * (p.x - o.x)) +
              o.z * (u.y * v.x - u.x * v.y) +
              p.z * (u.x * v.y - u.y * v.x)) / \
            (u.z * (v.x * w.y - v.y * w.x) +
             u.y * (v.z * w.x - v.x * w.z) +
             u.x * (v.y * w.z - v.z * w.y))

        return d1 >= 0 and d1 <= 1 and d2 >= 0 and d2 <= 1 and \
            d3 >= 0 and d3 <= 1


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

    def __repr__(self):
        return "Cuboid(pmin={!r}, pmax={!r})".format(
            self.pmin.tolist(), self.pmax.tolist())

    @property
    def pmin(self):
        return self._pmin

    @property
    def pmax(self):
        return self._pmax

    @property
    def xmin(self):
        return self.pmin.x

    @property
    def ymin(self):
        return self.pmin.y

    @property
    def zmin(self):
        return self.pmin.z

    @property
    def xmax(self):
        return self.pmax.x

    @property
    def ymax(self):
        return self.pmax.y

    @property
    def zmax(self):
        return self.pmax.z

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

    def contains_point(self, point):
        """Check if point is contained within volume of `Cuboid`."""
        p = Point(point)

        return (p.x >= self.xmin) and (p.x <= self.xmax) and \
            (p.y >= self.ymin) and (p.y <= self.ymax) and \
            (p.z >= self.zmin) and (p.z <= self.zmax)


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

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)

        self._center = center

        if a is None:
            a = 1.0
        self._a = a

    def __repr__(self):
        return("Cube(center={!r}, a={!r})".format(self.center.tolist(),
                                                  self.a))

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

    def contains_point(self, point):
        p = Point(point)
        c = self.center
        a = self.a
        xmin = c.x - a / 2
        ymin = c.y - a / 2
        zmin = c.z - a / 2
        xmax = c.x + a / 2
        ymax = c.y + a / 2
        zmax = c.z + a / 2
        return (p.x >= xmin) and (p.x <= xmax) and \
            (p.y >= ymin) and (p.y <= ymax) and \
            (p.z >= zmin) and (p.z <= zmax)


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

    def __repr__(self):
        return("Ellipsoid(center={!r}, a={!r}, b={!r}, c={!r})".format(
            self.center, self.a, self.b, self.c))

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

    def contains_point(self, point):
        """Check if point is contained within volume of :class:`Ellipsoid`."""
        p = Point(point)
        c = self.center
        a, b, c = self.a, self.b, self.c

        return (p.x - c.x)**2 / a**2 + (p.y - c.y)**2 / b**2 + \
            (p.z - c.z)**2 / c**2 <= 1.0


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

    def __repr__(self):
        return("Spheroid(center={!r}, a={!r}, c={!r})".format(
            self.center, self.a, self.c))

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

    def contains_point(self, point=None):
        """Check if point is contained within volume of :class:`Spheroid`."""
        x, y, z = point

        h, k, l = self.center
        a, c = self.a, self.c

        return ((x - h)**2 + (y - k)**2) / a**2 + (z - l)**2 / c**2 <= 1.0


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

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)

        self._center = center

        if r is None:
            r = 1.0

        self._r = r

    def __repr__(self):
        return("Sphere(center={!r}, r={!r})".format(self.center, self.r))

    @property
    def r(self):
        return self._r

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

    def contains_point(self, point=None):
        x, y, z = point

        h, k, l = self.center
        r = self.r

        return (x - h)**2 + (y - k)**2 + (z - l)**2 <= r**2
