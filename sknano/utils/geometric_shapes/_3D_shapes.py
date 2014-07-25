# -*- coding: utf-8 -*-
"""
=======================================================================
3D geometric shapes (:mod:`sknano.utils.geometric_shapes._3D_shapes`)
=======================================================================

.. currentmodule:: sknano.utils.geometric_shapes._3D_shapes

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

import numpy as np

from sknano.core import Point

__all__ = ['Cube', 'Cuboid', 'Ellipsoid', 'Spheroid', 'Sphere',
           'Polyhedron', 'Hexahedron', 'Parallelepiped', 'Rhombohedron']


class Parallelepiped(object):
    """Abstract base class for defining common properties of
    geometric shapes with 6 parallel sides.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    xmin, ymin, zmin : float
    xmax, ymax, zmax : float
    pmin, pmax : sequence or :class:`~sknano.core.Point`, optional

    """
    __metaclass__ = ABCMeta

    def __init__(self, xmin=None, ymin=None, zmin=None,
                 xmax=None, ymax=None, zmax=None, pmin=None, pmax=None):

        if pmin is None:
            pmin = Point(x=xmin, y=ymin, z=zmin)
        elif isinstance(pmin, (tuple, list, np.ndarray)):
            x, y, z = pmin
            pmin = Point(x=x, y=y, z=z)

        self._pmin = pmin
        self._xmin, self._ymin, self._zmin = self._pmin

        if pmax is None:
            pmax = Point(x=xmax, y=ymax, z=zmax)
        elif isinstance(pmax, (tuple, list, np.ndarray)):
            x, y, z = pmax
            pmax = Point(x=x, y=y, z=z)

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

    @abstractmethod
    def contains_point(self):
        """Check if point is contained within volume of parallelpiped."""
        return NotImplementedError('Subclasses of `Parallelepiped` must '
                                   'implement a `contains_point` method.')


class Cuboid(Parallelepiped):
    """Abstract data structure representing a cuboid.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    xmin, ymin, zmin : float
    xmax, ymax, zmax : float
    pmin, pmax : sequence or :class:`~sknano.core.Point`, optional

    """
    def __init__(self, xmin=None, ymin=None, zmin=None,
                 xmax=None, ymax=None, zmax=None, pmin=None, pmax=None):

        super(Cuboid, self).__init__(xmin=xmin, ymin=ymin, zmin=zmin,
                                     xmax=xmax, ymax=ymax, zmax=zmax,
                                     pmin=pmin, pmax=pmax)

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

    def contains_point(self, point=None):
        """Check if point is contained within volume of cuboid."""
        x, y, z = point

        return((x > self._xmin) and (x < self._xmax) and
               (y > self._ymin) and (y < self._ymax) and
               (z > self._zmin) and (z < self._zmax))


class Cube(Cuboid):
    """Abstract data structure representing a cube.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    a : float, optional
        length of edge
    center : sequence or :class:`~sknano.core.Point`
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the cube center.

    """
    def __init__(self, a=None, center=None):

        if a is None:
            raise TypeError('Please specify the edge length parameter `a`')

        self._a = a

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(x=center[0], y=center[1], z=center[2])

        self._center = center

        h, k, l = self._center
        super_kwargs = \
            {'xmin': h - a / 2, 'ymin': k - a / 2, 'zmin': l - a / 2,
             'xmax': h + a / 2, 'ymax': k + a / 2, 'zmax': l + a / 2}

        super(Cube, self).__init__(**super_kwargs)

    def __repr__(self):
        return("Cube(a={!r}, center={!r})".format(self.a, self.center))

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        raise AttributeError("{!r} object has no attribute {!r}".format(
            self.__class__.__name__, 'b'))

    @property
    def c(self):
        raise AttributeError("{!r} object has no attribute {!r}".format(
            self.__class__.__name__, 'c'))

    @property
    def center(self):
        return self._center


class Ellipsoid(object):
    """Abstract data structure representing an ellipsoid.

    .. versionadded:: 0.2.26

    The general ellipsoid is a quadratic surface with is given in
    Cartesian coordinates by:

    .. math::

       \\frac{x^2}{a^2} + \\frac{y^2}{b^2} + \\frac{z^2}{c^2} = 1

    Parameters
    ----------
    a, b, c : float, optional
        Semi-principal axes :math:`a, b, c` of axis-aligned
        :class:`Ellipsoid`
    center : sequence or :class:`~sknano.core.Point`
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Ellipsoid` center.

    """
    def __init__(self, a=None, b=None, c=None, center=None):

        if a is None or b is None or c is None:
            raise TypeError("Parameters 'a', 'b', 'c' must be specified.")

        self._a, self._b, self._c = a, b, c

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            h, k, l = center
            center = Point(x=h, y=k, z=l)

        self._center = center

    def __repr__(self):
        return("Ellipsoid(a={!r}, b={!r}, c={!r}, center={!r})".format(
            self.a, self.b, self.c, self.center))

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

    def contains_point(self, point=None):
        """Check if point is contained within volume of :class:`Ellipsoid`."""
        x, y, z = point

        h, k, l = self.center
        a, b, c = self.a, self.b, self.c

        return((x - h)**2 / a**2 +
               (y - k)**2 / b**2 +
               (z - l)**2 / c**2 < 1.0)


class Spheroid(object):
    """Abstract data structure representing a spheroid.

    .. versionadded:: 0.2.26

    The general spheroid is a quadratic surface with is given in
    Cartesian coordinates by:

    .. math::

       \\frac{x^2 + y^2}{a^2} + \\frac{z^2}{c^2} = 1

    Parameters
    ----------
    a, c : float, optional
        Semi-axes :math:`a, c` of axis-aligned
        :class:`Spheroid` with symmetry axis along :math:`z` axis.
    center : sequence or :class:`~sknano.core.Point`
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Spheroid` center.

    """
    def __init__(self, a=None, c=None, center=None):

        if a is None or c is None:
            raise TypeError("Parameters 'a', 'c' must be specified.")

        self._a, self._c = a, c

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            h, k, l = center
            center = Point(x=h, y=k, z=l)

        self._center = center

    def __repr__(self):
        return("Spheroid(a={!r}, c={!r}, center={!r})".format(
            self.a, self.c, self.center))

    @property
    def center(self):
        return self._center

    @property
    def a(self):
        return self._a

    @property
    def c(self):
        return self._c

    def contains_point(self, point=None):
        """Check if point is contained within volume of :class:`Ellipsoid`."""
        x, y, z = point

        h, k, l = self.center
        a, c = self.a, self.c

        return(((x - h)**2 + (y - k)**2) / a**2 +
               (z - l)**2 / c**2 < 1.0)


class Sphere(object):
    """Abstract data structure representing a sphere.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    r : float, optional
        Sphere radius :math:`r`
    center : sequence or :class:`~sknano.core.Point`
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Sphere` center.

    """
    def __init__(self, r=1.0, center=None):

        if r is None:
            raise TypeError('Please specify the sphere radius `r`')

        self._r = r

        if center is None:
            center = Point()
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = Point(x=center[0], y=center[1], z=center[2])

        self._center = center

    def __repr__(self):
        return("Sphere(r={!r}, center={!r})".format(self.r, self.center))

    @property
    def r(self):
        return self._r

    @property
    def center(self):
        return self._center

    def contains_point(self, point=None):
        x, y, z = point

        h, k, l = self.center
        r = self.r

        return((x - h)**2 + (y - k)**2 + (z - l)**2 < r**2)


class Polyhedron(object):
    pass


class Hexahedron(object):
    pass


class Rhombohedron(object):
    pass
