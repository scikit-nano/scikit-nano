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

from sknano.core import Point
from ._base import GeometricRegion

__all__ = ['Geometric3DRegion', 'Cube', 'Cuboid', 'Ellipsoid', 'Spheroid',
           'Sphere', 'Parallelepiped']


class Geometric3DRegion(GeometricRegion):
    """Abstract base class for representing 3D geometric regions."""
    __metaclass__ = ABCMeta

    @abstractproperty
    def volume(self):
        """Volume of 3D geometric region."""
        raise NotImplementedError


class Parallelepiped(Geometric3DRegion):
    """Abstract representation of parallelepiped.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    xmin, ymin, zmin : float
    xmax, ymax, zmax : float
    pmin, pmax : sequence, optional

    """
    def __init__(self, xmin=None, ymin=None, zmin=None,
                 xmax=None, ymax=None, zmax=None, pmin=None, pmax=None):

        if pmin is None:
            pmin = Point([xmin, ymin, zmin])
        elif isinstance(pmin, (tuple, list, np.ndarray)):
            pmin = Point(pmin)

        self._pmin = pmin
        self._xmin, self._ymin, self._zmin = self._pmin

        if pmax is None:
            pmax = Point([xmax, ymax, zmax])
        elif isinstance(pmax, (tuple, list, np.ndarray)):
            pmax = Point(pmax)

        self._pmax = pmax
        self._xmax, self._ymax, self._zmax = self._pmax

    @property
    def xmin(self):
        return self._xmin

    @property
    def ymin(self):
        return self._ymin

    @property
    def zmin(self):
        return self._zmin

    @property
    def xmax(self):
        return self._xmax

    @property
    def ymax(self):
        return self._ymax

    @property
    def zmax(self):
        return self._zmax

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
        l = (self.zmax + self.zmin) / 2
        return Point(x=h, y=k, z=l)

    @property
    def centroid(self):
        pass

    @property
    def volume(self):
        pass

    def contains_point(self):
        """Check if point is contained within volume of parallelpiped."""
        pass


class Cuboid(Geometric3DRegion):
    """Abstract data structure representing a cuboid.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    xmin, ymin, zmin : float
    xmax, ymax, zmax : float
    pmin, pmax : sequence, optional

    """
    def __init__(self, xmin=None, ymin=None, zmin=None,
                 xmax=None, ymax=None, zmax=None, pmin=None, pmax=None):

        pass

    def __repr__(self):
        return("Cuboid(xmin={!r}, ymin={!r}, zmin={!r}, "
               "xmax={!r}, ymax={!r}, zmax={!r})".format(
                   self._xmin, self._ymin, self._zmin,
                   self._xmax, self._ymax, self._zmax))

    @property
    def a(self):
        return self._xmax - self._xmin

    @property
    def b(self):
        return self._ymax - self._ymin

    @property
    def c(self):
        return self._zmax - self._zmin

    @property
    def centroid(self):
        pass

    @property
    def volume(self):
        pass

    def contains_point(self, point=None):
        """Check if point is contained within volume of cuboid."""
        x, y, z = point

        return((x > self._xmin) and (x < self._xmax) and
               (y > self._ymin) and (y < self._ymax) and
               (z > self._zmin) and (z < self._zmax))


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
    def __init__(self, a=None, center=None):

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)

        self._center = center

        if a is None:
            a = 1.0

        self._a = a

        h, k, l = self._center

        bounds = \
            {'xmin': h - a / 2, 'ymin': k - a / 2, 'zmin': l - a / 2,
             'xmax': h + a / 2, 'ymax': k + a / 2, 'zmax': l + a / 2}

    def __repr__(self):
        return("Cube(center={!r}, a={!r})".format(self.center, self.a))

    @property
    def a(self):
        return self._a

    @property
    def center(self):
        return self._center

    @property
    def centroid(self):
        pass

    @property
    def volume(self):
        pass

    def contains_point(self, point=None):
        pass


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

        if a is None:
            a = 0.5
        if b is None:
            b = 0.5
        if c is None:
            c = 0.5

        self._a, self._b, self._c = a, b, c

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(center)

        self._center = center

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
        pass

    @property
    def volume(self):
        pass

    def contains_point(self, point=None):
        """Check if point is contained within volume of :class:`Ellipsoid`."""
        x, y, z = point

        h, k, l = self.center
        a, b, c = self.a, self.b, self.c

        return((x - h)**2 / a**2 +
               (y - k)**2 / b**2 +
               (z - l)**2 / c**2 < 1.0)


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
        pass

    @property
    def volume(self):
        pass

    def contains_point(self, point=None):
        """Check if point is contained within volume of :class:`Ellipsoid`."""
        x, y, z = point

        h, k, l = self.center
        a, c = self.a, self.c

        return(((x - h)**2 + (y - k)**2) / a**2 +
               (z - l)**2 / c**2 < 1.0)


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
    def __init__(self, center=None, r=1.0):

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
        pass

    @property
    def volume(self):
        pass

    def contains_point(self, point=None):
        x, y, z = point

        h, k, l = self.center
        r = self.r

        return((x - h)**2 + (y - k)**2 + (z - l)**2 < r**2)
