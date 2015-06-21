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
           'Ellipse', 'Circle', 'Triangle']


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
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._o = Point(value)

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._u = Vector(value, p0=self.o)

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
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

        cx = 0.5 * (2 * ox + ux + vx)
        cy = 0.5 * (2 * oy + uy + vy)

        return Point([cx, cy])

    def contains(self, point):
        """Check if `point` is within region."""
        px, py = Point(point)

        ox, oy = self.o
        ux, uy = self.u
        vx, vy = self.v

        q1 = ((py - oy) * vx + (ox - px) * vy) / (uy * vx - ux * vy)
        q2 = ((py - oy) * ux + (ox - px) * uy) / (ux * vy - uy * vx)

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
    def __init__(self, pmin=None, pmax=None, xmin=0, ymin=0, xmax=1, ymax=1):

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
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._pmin = Point(value)

    @property
    def pmax(self):
        return self._pmax

    @pmax.setter
    def pmax(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
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
        cx = (self.xmax + self.xmin) / 2
        cy = (self.ymax + self.ymin) / 2
        return Point([cx, cy])

    def contains(self, point):
        """Check if `point` is within region."""
        px, py = Point(point)
        xmin = self.xmin
        xmax = self.xmax
        ymin = self.ymin
        ymax = self.ymax

        return (px >= xmin) and (px <= xmax) and (py >= ymin) and (py <= ymax)

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
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
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
        return self.a ** 2

    @property
    def centroid(self):
        return self.center

    def contains(self, point):
        px, py = Point(point)
        cx, cy = self.center
        a = self.a
        xmin = cx - a / 2
        xmax = cx + a / 2
        ymin = cy - a / 2
        ymax = cy + a / 2
        return (px >= xmin) and (px <= xmax) and (py >= ymin) and (py <= ymax)

    def todict(self):
        return dict(center=self.center, a=self.a)


class Triangle(Geometric2DRegion):
    """`Geometric2DRegion` for a triangle.

    .. versionadded:: 0.3.10

    Parameters
    ----------
    p1, p2, p3 : sequence, optional
        2-tuples or :class:`~sknano.core.Point` class instances
        specifying the `Triangle` corner points :math:`p_1=(x_1,y_1)`,
        :math:`p_2=(x_2,y_2)`, and :math:`p_3=(x_3, y_3)`.
    """
    def __init__(self, p1=None, p2=None, p3=None):
        super().__init__()

        if p1 is None:
            p1 = [0, 0]
        if p2 is None:
            p2 = [0, 1]
        if p3 is None:
            p3 = [1, 0]

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.points.extend([self.p1, self.p2, self.p3])
        self.fmtstr = "p1={p1!r}, p2={p2!r}, p3={p3!r}"

    @property
    def p1(self):
        return self._p1

    @p1.setter
    def p1(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._p1 = Point(value)

    @property
    def p2(self):
        return self._p2

    @p2.setter
    def p2(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._p2 = Point(value)

    @property
    def p3(self):
        return self._p3

    @p3.setter
    def p3(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._p3 = Point(value)

    @property
    def centroid(self):
        x1, y1 = self.p1
        x2, y2 = self.p2
        x3, y3 = self.p3

        cx = (x1 + x2 + x3) / 3
        cy = (y1 + y2 + y3) / 3
        return Point([cx, cy])

    @property
    def area(self):
        x1, y1 = self.p1
        x2, y2 = self.p2
        x3, y3 = self.p3

        return np.abs(-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 +
                      x2 * y3) / 2

    def contains(self, point):
        px, py = Point(point)
        x1, y1 = self.p1
        x2, y2 = self.p2
        x3, y3 = self.p3

        q1 = ((x1 - x3) * py + (x3 - px) * y1 + (px - x1) * y3) / \
            ((y1 - y2) * x3 + (y2 - y3) * x1 + (y3 - y1) * x2)

        q2 = ((x2 - x1) * py + (px - x2) * y1 + (x1 - px) * y2) / \
            ((y1 - y2) * x3 + (y2 - y3) * x1 + (y3 - y1) * x2)

        q3 = ((x2 - x3) * py + (x3 - px) * y2 + (px - x2) * y3) / \
            ((y1 - y2) * x3 + (y2 - y3) * x1 + (y3 - y1) * x2)

        return q1 >= 0 and q2 >= 0 and q3 <= 0

    def todict(self):
        return dict(p1=self.p1, p2=self.p2, p3=self.p3)


class Ellipse(Geometric2DRegion):
    """Abstract data structure representing an ellipse.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
        Center of axis-aligned ellipse with semi-axes lengths :math:`r_x, r_y`
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
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
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
    def area(self):
        return np.pi * self.rx * self.ry

    @property
    def centroid(self):
        return self.center

    def contains(self, point):
        px, py = Point(point)
        cx, cy = self.center
        rx, ry = self.rx, self.ry

        q1 = (px - cx) ** 2 / rx ** 2
        q2 = (py - cy) ** 2 / ry ** 2

        return q1 + q2 <= 1.0

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
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
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
        return np.pi * self.r ** 2

    @property
    def centroid(self):
        return self.center

    def contains(self, point):
        x, y = Point(point)
        h, k = self.center
        r = self.r

        return (x - h) ** 2 + (y - k) ** 2 <= r ** 2

    def todict(self):
        return dict(center=self.center, r=self.r)
