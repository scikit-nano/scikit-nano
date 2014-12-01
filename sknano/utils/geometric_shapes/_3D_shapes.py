# -*- coding: utf-8 -*-
"""
=======================================================================
3D geometric shapes (:mod:`sknano.utils.geometric_shapes._3D_shapes`)
=======================================================================

.. currentmodule:: sknano.utils.geometric_shapes._3D_shapes

"""
from __future__ import absolute_import, division, print_function
import six
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractproperty

import numbers
import numpy as np

from sknano.core.math import Point, Vector, vector as vec
from ._base import GeometricRegion, GeometricTransformsMixin

__all__ = ['Geometric3DRegion', 'Parallelepiped', 'Cuboid', 'Cube',
           'Ellipsoid', 'Spheroid', 'Sphere']


class Geometric3DRegion(six.with_metaclass(ABCMeta, GeometricTransformsMixin,
                                           GeometricRegion)):
    """Abstract base class for representing 3D geometric regions."""

    @abstractproperty
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
        super(Parallelepiped, self).__init__()

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

    def update_region_limits(self):
        self.limits['x']['min'] = self.xmin
        self.limits['x']['max'] = self.xmax
        self.limits['y']['min'] = self.ymin
        self.limits['y']['max'] = self.ymax
        self.limits['z']['min'] = self.zmin
        self.limits['z']['max'] = self.zmax


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
        super(Cuboid, self).__init__()

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

    def __repr__(self):
        return "Cuboid(pmin={!r}, pmax={!r})".format(
            self.pmin.tolist(), self.pmax.tolist())

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

    def contains_point(self, point):
        """Check if point is contained within volume of `Cuboid`."""
        p = Point(point)

        return (p.x >= self.xmin) and (p.x <= self.xmax) and \
            (p.y >= self.ymin) and (p.y <= self.ymax) and \
            (p.z >= self.zmin) and (p.z <= self.zmax)

    def update_region_limits(self):
        self.limits['x']['min'] = self.xmin
        self.limits['x']['max'] = self.xmax
        self.limits['y']['min'] = self.ymin
        self.limits['y']['max'] = self.ymax
        self.limits['z']['min'] = self.zmin
        self.limits['z']['max'] = self.zmax


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
        super(Cube, self).__init__()

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if a is None:
            a = 1.0
        self._a = a

        self.points.append(self.center)

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

    def update_region_limits(self):
        self.limits['x']['min'] = self.xmin
        self.limits['x']['max'] = self.xmax
        self.limits['y']['min'] = self.ymin
        self.limits['y']['max'] = self.ymax
        self.limits['z']['min'] = self.zmin
        self.limits['z']['max'] = self.zmax


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
        super(Ellipsoid, self).__init__()

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

    def update_region_limits(self):
        c = self.center
        r = self.r
        self.limits['x']['min'] = c.x - r
        self.limits['x']['max'] = c.x + r
        self.limits['y']['min'] = c.y - r
        self.limits['y']['max'] = c.y + r
        self.limits['z']['min'] = c.z - r
        self.limits['z']['max'] = c.z + r


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
        super(Spheroid, self).__init__()

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

    def update_region_limits(self):
        c = self.center
        r = self.r
        self.limits['x']['min'] = c.x - r
        self.limits['x']['max'] = c.x + r
        self.limits['y']['min'] = c.y - r
        self.limits['y']['max'] = c.y + r
        self.limits['z']['min'] = c.z - r
        self.limits['z']['max'] = c.z + r


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
        super(Sphere, self).__init__()

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if r is None:
            r = 1.0
        self._r = r

        self.points.append(self.center)

    def __repr__(self):
        return("Sphere(center={!r}, r={!r})".format(self.center, self.r))

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

    def contains_point(self, point):
        p = point
        c = self.center
        r = self.r

        return (p.x - c.x)**2 + (p.y - c.y)**2 + (p.z - c.z)**2 <= r**2

    def update_region_limits(self):
        c = self.center
        r = self.r
        self.limits['x']['min'] = c.x - r
        self.limits['x']['max'] = c.x + r
        self.limits['y']['min'] = c.y - r
        self.limits['y']['max'] = c.y + r
        self.limits['z']['min'] = c.z - r
        self.limits['z']['max'] = c.z + r
