# -*- coding: utf-8 -*-
"""
======================================================================
2D geometric shapes (:mod:`sknano.utils.geometric_shapes._2D_shapes`)
======================================================================

.. currentmodule:: sknano.utils.geometric_shapes._2D_shapes

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractproperty

#import numbers
import numpy as np

from sknano.core import Point, Vector
from ._base import GeometricRegion

__all__ = ['Geometric2DRegion', 'Circle', 'Ellipse', 'Parallelogram',
           'Rhombus', 'Rhomboid', 'Rectangle', 'Square']


class Geometric2DRegion(GeometricRegion):
    """Abstract base class for representing 2D geometric regions."""
    __metaclass__ = ABCMeta

    @abstractproperty
    def area(self):
        """Area of 2D geometric region."""
        raise NotImplementedError


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

        if center is None or not isinstance(center, (tuple, list, np.ndarray)):
            center = Point(nd=2)
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if r is None:
            r = 1.0
        self._r = r

    def __repr__(self):
        return("Circle(center={!r}, r={!r})".format(self.center, self.r))

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
    def area(self):
        pass

    @property
    def centroid(self):
        pass

    def contains_point(self, point=None):
        h, k = self.center
        r = self.r

        x, y = point

        return (x - h)**2 + (y - k)**2 <= r**2


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

        if center is None or not isinstance(center, (tuple, list, np.ndarray)):
            center = Point(nd=2)
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)

        self._center = center
        self._rx = rx
        self._ry = ry

        self._a = rx
        self._b = ry
        if rx < ry:
            self._a = ry
            self._b = rx

    def __repr__(self):
        return("Ellipse(center={!r}, rx={!r}, ry={!r})".format(
            self.center, self.rx, self.ry))

    @property
    def center(self):
        return self._center

    @property
    def rx(self):
        return self._rx

    @property
    def ry(self):
        return self._ry

    @property
    def a(self):
        """Semi-major axis length"""
        return self._a

    @property
    def b(self):
        """Semi-minor axis length"""
        return self._b

    @property
    def area(self):
        pass

    @property
    def centroid(self):
        pass

    def contains_point(self, point=None):
        h, k = self.center
        rx, ry = self.rx, self.ry

        x, y = point

        return (x - h)**2 / rx**2 + (y - k)**2 / ry**2 <= 1.0


class Parallelogram(Geometric2DRegion):
    """Abstract object representation of a parallelogram.

    Represents a parallelogram with origin :math:`p` and directions
    :math:`v_1` and :math:`v_2`.

    Parameters
    ----------
    p : sequence, optional
        parallelogram origin
    v1, v2 : array_like, optional
        parallelogram direction vectors stemming from origin `p`.

    """
    def __init__(self, p=None, v1=None, v2=None):

        if p is None:
            p = Point(nd=2)
        elif isinstance(p, (tuple, list, np.ndarray)):
            p = Point(p)
        self._p = p

        if v1 is None or v2 is None:
            v1 = Vector([1., 0.])
            v2 = Vector([1., 1.])


class Rhombus(object):
    pass


class Rhomboid(object):
    pass


class Rectangle(Geometric2DRegion):
    """Abstract data structure representing a rectangle.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    xmin, ymin : float
    xmax, ymax : float
    pmin, pmax : sequence, optional

    """
    def __init__(self, xmin=None, ymin=None, xmax=None, ymax=None,
                 pmin=None, pmax=None):

        if pmin is None:
            pmin = Point([xmin, ymin])
        elif isinstance(pmin, (tuple, list, np.ndarray)):
            pmin = Point(pmin)

        self._pmin = pmin
        self._xmin, self._ymin = self._pmin

        if pmax is None:
            pmax = Point([xmax, ymax])
        elif isinstance(pmax, (tuple, list, np.ndarray)):
            pmax = Point(pmax)

        self._pmax = pmax
        self._xmax, self._ymax = self._pmax

    def __repr__(self):
        return("Rectangle(xmin={!r}, ymin={!r}, xmax={!r}, ymax={!r})".format(
            self._xmin, self._ymin, self._xmax, self._ymax))

    @property
    def xmin(self):
        return self._xmin

    @property
    def ymin(self):
        return self._ymin

    @property
    def xmax(self):
        return self._xmax

    @property
    def ymax(self):
        return self._ymax

    @property
    def pmin(self):
        return self._pmin

    @property
    def pmax(self):
        return self._pmax

    @property
    def center(self):
        h = (self.xmax + self.xmin) / 2
        k = (self.ymax + self.ymin) / 2
        return Point([h, k])

    @property
    def a(self):
        return self._xmax - self._xmin

    @property
    def b(self):
        return self._ymax - self._ymin

    @property
    def area(self):
        pass

    @property
    def centroid(self):
        pass

    def contains_point(self, point=None):
        """Check if point is contained within volume of cuboid."""
        x, y = point

        return (x >= self._xmin) and (x <= self._xmax) and \
            (y >= self._ymin) and (y <= self._ymax)


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

        if center is None:
            center = Point(nd=2)
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)
        self._center = center

        if a is None:
            a = 1.0
        self._a = a

    def __repr__(self):
        return("Square(center={!r}, a={!r})".format(self.center, self.a))

    @property
    def center(self):
        return self._center

    @property
    def a(self):
        return self._a

    @property
    def area(self):
        pass

    @property
    def centroid(self):
        pass

    def contains_point(self, point=None):
        pass
