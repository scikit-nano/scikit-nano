# -*- coding: utf-8 -*-
"""
======================================================================
2D geometric shapes (:mod:`sknano.tools.geometric_shapes._2D_shapes`)
======================================================================

.. currentmodule:: sknano.tools.geometric_shapes._2D_shapes

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

import numpy as np

from sknano.tools import Point2D

__all__ = ['Circle', 'Ellipse', 'Polygon', 'Parallelogram', 'Rhombus',
           'Rhomboid', 'Rectangle', 'Square']


class Circle(object):
    """Abstract data structure representing a circle.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    r : float
        Circle radius.
    center : sequence or :class:`~sknano.tools.Point2D`
        Center of circle

    """
    def __init__(self, r=1.0, center=None):

        if center is None:
            center = Point2D()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point2D(x=center[0], y=center[1])

        self._r = r
        self._center = center

    def __repr__(self):
        return("Circle(r={!r}, center={!r})".format(self.r, self.center))

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value

    @property
    def center(self):
        return self._center

    def contains_point(self, point=None):
        h, k = self.center
        r = self.r

        x, y = point

        return((x - h)**2 + (y - k)**2 < r**2)


class Ellipse(object):
    """Abstract data structure representing an ellipse.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    center : sequence or :class:`~sknano.tools.Point2D`
        Center of axis-aligned ellipse with semi-axes :math:`r_x, r_y`
    rx, ry : float
        Lengths of semi-axes :math:`r_x, r_y`

    """
    def __init__(self, center=None, rx=1, ry=1):

        if center is None:
            center = Point2D()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point2D(x=center[0], y=center[1])

        self._center = center
        self._rx = rx
        self._ry = ry

        self._a = rx if rx >= ry else ry
        self._b = ry if rx >= ry else rx

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

    def contains_point(self, point=None):
        h, k = self.center
        rx, ry = self.rx, self.ry

        x, y = point

        return((x - h)**2 / rx**2 + (y - k)**2 / ry**2 < 1.0)


class Polygon(object):
    pass


class Parallelogram(object):
    """Abstract base class for defining common properties of geometric shapes
    that are parallelograms or a subset of.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    xmin, ymin : float
    xmax, ymax : float
    pmin, pmax : sequence or :class:`~sknano.tools.Point2D`, optional

    """
    __metaclass__ = ABCMeta

    def __init__(self, xmin=None, ymin=None, xmax=None, ymax=None,
                 pmin=None, pmax=None):

        if pmin is None:
            pmin = Point2D(x=xmin, y=ymin)
        elif isinstance(pmin, (tuple, list, np.ndarray)):
            x, y = pmin
            pmin = Point2D(x=x, y=y)

        self._pmin = pmin
        self._xmin, self._ymin = self._pmin

        if pmax is None:
            pmax = Point2D(x=xmax, y=ymax)
        elif isinstance(pmax, (tuple, list, np.ndarray)):
            x, y = pmax
            pmax = Point2D(x=x, y=y)

        self._pmax = pmax
        self._xmax, self._ymax = self._pmax

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
        return Point2D(x=h, y=k)

    @abstractmethod
    def contains_point(self):
        """Check if point is contained within volume of parallelogram."""
        return NotImplementedError('Subclasses of `Parallelogram` must '
                                   'implement a `contains_point` method.')


class Rhombus(object):
    pass


class Rhomboid(object):
    pass


class Rectangle(Parallelogram):
    """Abstract data structure representing a rectangle.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    xmin, ymin : float
    xmax, ymax : float
    pmin, pmax : sequence or :class:`~sknano.tools.Point2D`, optional

    """
    def __init__(self, xmin=None, ymin=None, xmax=None, ymax=None,
                 pmin=None, pmax=None):

        super(Rectangle, self).__init__(xmin=xmin, ymin=ymin,
                                        xmax=xmax, ymax=ymax,
                                        pmin=pmin, pmax=pmax)

    def __repr__(self):
        return("Rectangle(xmin={!r}, ymin={!r}, xmax={!r}, ymax={!r})".format(
            self._xmin, self._ymin, self._xmax, self._ymax))

    @property
    def a(self):
        return self._xmax - self._xmin

    @property
    def b(self):
        return self._ymax - self._ymin

    def contains_point(self, point=None):
        """Check if point is contained within volume of cuboid."""
        x, y = point

        return((x > self._xmin) and (x < self._xmax) and
               (y > self._ymin) and (y < self._ymax))


class Square(Rectangle):
    """Abstract data structure representing a square.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    a : float, optional
        length of side
    center : sequence or :class:`~sknano.tools.Point2D`, optional

    """
    def __init__(self, a=None, center=None):

        if a is None:
            raise TypeError('Please specify the edge length parameter `a`')

        self._a = a

        if center is None:
            center = Point2D()
        elif isinstance(center, (tuple, list, np.ndarray)):
            h, k = center
            center = Point2D(x=h, y=k)

        self._center = center

        h, k = self._center
        super_kwargs = {'xmin': h - a / 2, 'ymin': k - a / 2,
                        'xmax': h + a / 2, 'ymax': k + a / 2}

        super(Square, self).__init__(**super_kwargs)

    def __repr__(self):
        return("Square(a={!r}, center={!r})".format(self.a, self.center))

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        raise AttributeError("{!r} object has no attribute {!r}".format(
            self.__class__.__name__, 'b'))

    @property
    def center(self):
        return self._center
