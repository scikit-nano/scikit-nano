# -*- coding: utf-8 -*-
"""
======================================================================
2D geometric shapes (:mod:`sknano.utils.geometric_shapes._2D_shapes`)
======================================================================

.. currentmodule:: sknano.utils.geometric_shapes._2D_shapes

"""
from __future__ import absolute_import, division, print_function
import six
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractproperty

import numpy as np

from sknano.core.math import Point, Vector
from ._base import GeometricRegion, GeometricTransformsMixin

__all__ = ['Geometric2DRegion', 'Parallelogram', 'Rectangle', 'Square',
           'Ellipse', 'Circle']


class Geometric2DRegion(six.with_metaclass(ABCMeta, GeometricTransformsMixin,
                                           GeometricRegion)):
    """Abstract base class for representing 2D geometric regions."""

    @abstractproperty
    def area(self):
        """Area of 2D geometric region."""
        raise NotImplementedError


class Parallelogram(Geometric2DRegion):
    """Abstract object representation of a parallelogram.

    Represents a parallelogram with origin :math:`o` and directions
    :math:`u` and :math:`v`.

    Parameters
    ----------
    o : array_like, optional
        parallelogram origin
    u, v : array_like, optional
        parallelogram direction vectors stemming from origin `o`.

    """
    def __init__(self, o=None, u=None, v=None):

        super(Parallelogram, self).__init__()

        if o is None:
            o = Point(nd=2)
        elif isinstance(o, (tuple, list, np.ndarray)):
            o = Point(o)
        self._o = o

        if u is None:
            u = [1., 0.]
        self._u = Vector(u, p0=self._o)

        if v is None:
            v = [1., 1.]
        self._v = Vector(v, p0=self._o)

        self.points.append(self.o)
        self.vectors.extend([self.u, self.v])

    def __repr__(self):
        return "Parallelogram(o={!r}, u={!r}, v={!r})".format(
            self.o.tolist(), self.u, self.v)

    @property
    def o(self):
        return self._o

    @property
    def u(self):
        return self._u

    @property
    def v(self):
        return self._v

    @property
    def center(self):
        return self.centroid

    @property
    def area(self):
        u = self.u
        v = self.v
        return np.abs(np.cross(u, v))

    @property
    def centroid(self):
        o = self.o
        u = self.u
        v = self.v

        xcom = 0.5 * (2 * o.x + u.x + v.x)
        ycom = 0.5 * (2 * o.y + u.y + v.y)

        return Point([xcom, ycom])

    def contains_point(self, point):
        """Check if point is contained within volume of cuboid."""
        p = Point(point)

        o = self.o
        u = self.u
        v = self.v

        d1 = ((p.y - o.y) * v.x + (o.x - p.x) * v.y) / (u.y * v.x - u.x * v.y)
        d2 = ((p.y - o.y) * u.x + (o.x - p.x) * u.y) / (u.x * v.y - u.y * v.x)

        return d1 >= 0 and d1 <= 1 and d2 >= 0 and d2 <= 1


class Rectangle(Geometric2DRegion):
    """Abstract data structure representing a rectangle.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    xmin, ymin : float
    xmax, ymax : float
    pmin, pmax : sequence, optional

    """
    def __init__(self, pmin=None, pmax=None, xmin=None, ymin=None,
                 xmax=None, ymax=None):

        super(Rectangle, self).__init__()

        if pmin is None:
            pmin = Point([xmin, ymin])
        elif isinstance(pmin, (tuple, list, np.ndarray)):
            pmin = Point(pmin)
        self._pmin = pmin

        if pmax is None:
            pmax = Point([xmax, ymax])
        elif isinstance(pmax, (tuple, list, np.ndarray)):
            pmax = Point(pmax)
        self._pmax = pmax

        self.points.append([self.pmin, self.pmax])

    def __repr__(self):
        return "Rectangle(pmin={!r}, pmax={!r})".format(
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
    def xmax(self):
        return self.pmax.x

    @property
    def ymin(self):
        return self.pmin.y

    @property
    def ymax(self):
        return self.pmax.y

    @property
    def center(self):
        return self.centroid

    @property
    def a(self):
        return self.xmax - self.xmin

    @property
    def b(self):
        return self.ymax - self.ymin

    @property
    def area(self):
        a = self.a
        b = self.b
        return a * b

    @property
    def centroid(self):
        h = (self.xmax + self.xmin) / 2
        k = (self.ymax + self.ymin) / 2
        return Point([h, k])

    def contains_point(self, point):
        """Check if point is contained within volume of cuboid."""
        p = Point(point)
        xmin = self.xmin
        xmax = self.xmax
        ymin = self.ymin
        ymax = self.ymax

        return (p.x >= xmin) and (p.x <= xmax) and \
            (p.y >= ymin) and (p.y <= ymax)


class Square(Geometric2DRegion):
    """Abstract data structure representing a square.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
    a : float, optional
        length of side

    """
    def __init__(self, center=None, a=None):

        super(Square, self).__init__()

        if center is None:
            center = Point(nd=2)
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if a is None:
            a = 1.0
        self._a = a

        self.points.append(self.center)

    def __repr__(self):
        return "Square(center={!r}, a={!r})".format(
            self.center.tolist(), self.a)

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        self._center = Point(value)

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def area(self):
        a = self.a
        return a**2

    @property
    def centroid(self):
        return self.center

    def contains_point(self, point):
        p = Point(point)
        c = self.center
        a = self.a
        xmin = c.x - a / 2
        ymin = c.y - a / 2
        xmax = c.x + a / 2
        ymax = c.y + a / 2
        return (p.x >= xmin) and (p.x <= xmax) and \
            (p.y >= ymin) and (p.y <= ymax)


class Ellipse(Geometric2DRegion):
    """Abstract data structure representing an ellipse.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
        Center of axis-aligned ellipse with semi-axes :math:`r_x, r_y`
    rx, ry : float
        Lengths of semi-axes :math:`r_x, r_y`

    """
    def __init__(self, center=None, rx=1, ry=1):
        super(Ellipse, self).__init__()

        if center is None or not isinstance(center, (tuple, list, np.ndarray)):
            center = Point(nd=2)
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)

        self._center = center
        self._rx = rx
        self._ry = ry

        self.points.append(self.center)

    def __repr__(self):
        return "Ellipse(center={!r}, rx={!r}, ry={!r})".format(
            self.center.tolist(), self.rx, self.ry)

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        self._center = Point(value)

    @property
    def rx(self):
        return self._rx

    @rx.setter
    def rx(self, value):
        self._rx = value

    @property
    def ry(self):
        return self._ry

    @ry.setter
    def ry(self, value):
        self._ry = value

    @property
    def a(self):
        """Semi-major axis length"""
        rx = self.rx
        ry = self.ry
        if rx < ry:
            return ry
        else:
            return rx

    @property
    def b(self):
        """Semi-minor axis length"""
        rx = self.rx
        ry = self.ry
        if rx < ry:
            return rx
        else:
            return ry

    @property
    def area(self):
        a = self.a
        b = self.b
        return np.pi * a * b

    @property
    def centroid(self):
        return self.center

    def contains_point(self, point):
        p = Point(point)
        c = self.center
        rx, ry = self.rx, self.ry

        return (p.x - c.x)**2 / rx**2 + (p.y - c.y)**2 / ry**2 <= 1.0


class Circle(Geometric2DRegion):
    """Abstract data structure representing a circle.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
        Center of circle
    r : float, optional
        Circle radius.

    """
    def __init__(self, center=None, r=1.0):

        super(Circle, self).__init__()

        if center is None or not isinstance(center, (tuple, list, np.ndarray)):
            center = Point(nd=2)
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if r is None:
            r = 1.0
        self._r = r

        self.points.append(self.center)

    def __repr__(self):
        return "Circle(center={!r}, r={!r})".format(
            self.center.tolist(), self.r)

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        self._center = Point(value)

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value

    @property
    def area(self):
        r = self.r
        return np.pi * r**2

    @property
    def centroid(self):
        return self.center

    def contains_point(self, point):
        p = Point(point)
        c = self.center
        r = self.r

        return (p.x - c.x)**2 + (p.y - c.y)**2 <= r**2
