# -*- coding: utf-8 -*-
"""
=======================================================================
2D geometric regions (:mod:`sknano.core.geometric_regions._2D_regions`)
=======================================================================

.. currentmodule:: sknano.core.geometric_regions._2D_regions

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

import numpy as np

from sknano.core.math import Point, Vector
from ._base import GeometricRegion, GeometricTransformsMixin

__all__ = ['Geometric2DRegion', 'Parallelogram', 'Rectangle', 'Square',
           'Ellipse', 'Circle']


class Geometric2DRegion(GeometricRegion, GeometricTransformsMixin,
                        metaclass=ABCMeta):
    """Abstract base class for representing 2D geometric regions."""

    @property
    @abstractmethod
    def area(self):
        """Area of 2D geometric region."""
        raise NotImplementedError

    @property
    def measure(self):
        """Measure of 2D geometric region."""
        return self.area


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

        super().__init__()

        if o is None:
            o = [0, 0]
        self.o = o

        if u is None:
            u = [1, 0]
        self.u = u

        if v is None:
            v = [1, 1]
        self.v = v

        self.points.append(self.o)
        self.vectors.extend([self.u, self.v])
        self.fmtstr = "o={o!r}, u={u!r}, v={v!r}"

    @property
    def o(self):
        return self._o

    @o.setter
    def o(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._o = Point(value)

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._u = Vector(value, p0=self.o)

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._v = Vector(value, p0=self.o)

    @property
    def area(self):
        u = self.u
        v = self.v
        return np.abs(np.cross(u, v))

    @property
    def centroid(self):
        ox, oy = self.o
        ux, uy = self.u
        vx, vy = self.v

        xcom = 0.5 * (2 * ox + ux + vx)
        ycom = 0.5 * (2 * oy + uy + vy)

        return Point([xcom, ycom])

    def contains(self, point):
        """Check if `point` is within region."""
        x, y = Point(point)

        ox, oy = self.o
        ux, uy = self.u
        vx, vy = self.v

        q1 = ((y - oy) * vx + (ox - x) * vy) / (uy * vx - ux * vy)
        q2 = ((y - oy) * ux + (ox - x) * uy) / (ux * vy - uy * vx)

        return q1 >= 0 and q1 <= 1 and q2 >= 0 and q2 <= 1

    def todict(self):
        return dict(o=self.o, u=self.u, v=self.v)


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

        super().__init__()

        if pmin is None:
            pmin = [xmin, ymin]
        self.pmin = pmin

        if pmax is None:
            pmax = [xmax, ymax]
        self.pmax = pmax

        self.points.append([self.pmin, self.pmax])

        self.fmtstr = "pmin={pmin!r}, pmax={pmax!r}"

    @property
    def pmin(self):
        return self._pmin

    @pmin.setter
    def pmin(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._pmin = Point(value)

    @property
    def pmax(self):
        return self._pmax

    @pmax.setter
    def pmax(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._pmax = Point(value)

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

    def contains(self, point):
        """Check if `point` is within region."""
        x, y = Point(point)
        xmin = self.xmin
        xmax = self.xmax
        ymin = self.ymin
        ymax = self.ymax

        return (x >= xmin) and (x <= xmax) and (y >= ymin) and (y <= ymax)

    def todict(self):
        return dict(pmin=self.pmin, pmax=self.pmax)


class Square(Geometric2DRegion):
    """Abstract data structure representing a square.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
    a : float, optional
        length of side

    """
    def __init__(self, center=None, a=1):

        super().__init__()

        if center is None:
            center = [0, 0]
        self.center = center

        self.a = a

        self.points.append(self.center)
        self.fmtstr = "center={center!r}, a={a:.2f}"

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
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

    def contains(self, point):
        x, y = Point(point)
        h, k = self.center
        a = self.a
        xmin = h - a / 2
        ymin = k - a / 2
        xmax = h + a / 2
        ymax = k + a / 2
        return (x >= xmin) and (x <= xmax) and (y >= ymin) and (y <= ymax)

    def todict(self):
        return dict(center=self.center, a=self.a)


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
        super().__init__()

        if center is None:
            center = [0, 0]
        self.center = center
        self.rx = rx
        self.ry = ry

        self.points.append(self.center)
        self.fmtstr = "center={center!r}, rx={rx:.3f}, ry={ry:.3f}"

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
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

    def contains(self, point):
        x, y = Point(point)
        h, k = self.center
        rx, ry = self.rx, self.ry

        return (x - h)**2 / rx**2 + (y - k)**2 / ry**2 <= 1.0

    def todict(self):
        return dict(center=self.center, rx=self.rx, ry=self.ry)


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

        super().__init__()

        if center is None:
            center = [0, 0]
        self.center = center
        self.r = r

        self.points.append(self.center)
        self.fmtstr = "center={center!r}, r={r:.3f}"

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
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

    def contains(self, point):
        x, y = Point(point)
        h, k = self.center
        r = self.r

        return (x - h)**2 + (y - k)**2 <= r**2

    def todict(self):
        return dict(center=self.center, r=self.r)
