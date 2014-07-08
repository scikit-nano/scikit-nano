# -*- coding: utf-8 -*-
"""
========================================================================
Classes to desribe 3D volume bounds (:mod:`sknano.nanogen._polyhedra`)
========================================================================

.. currentmodule:: sknano.structure_io.atoms._atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#from collections import OrderedDict

#import numpy as np

from ..tools import Point

__all__ = ['Cube', 'Cuboid']

_xyz = ('x', 'y', 'z')


class Cuboid(object):
    """Class defining a 'cuboid'.

    Parameters
    ----------
    xmin, ymin, zmin : float
    xmax, ymax, zmax : float

    """
    def __init__(self, xmin=None, ymin=None, zmin=None,
                 xmax=None, ymax=None, zmax=None):
        self._xmin = xmin
        self._ymin = ymin
        self._zmin = zmin
        self._xmax = xmax
        self._ymax = ymax
        self._zmax = zmax

    def __repr__(self):
        return("Cuboid(xmin={!r}, ymin={!r}, zmin={!r}, "
               "xmax={!r}, ymax={!r}, zmax={!r})".format(
                   self.xmin, self.ymin, self.zmin,
                   self.xmax, self.ymax, self.zmax))

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

    def contains_point(self, point=None):
        """Check if point is contained within volume of cuboid."""
        pt = None
        try:
            pt = Point(x=point[0], y=point[1], z=point[2])
        except (IndexError, TypeError):
            pt = Point

        try:
            if (pt.x > self.xmin) and (pt.x < self.xmax) and \
                    (pt.y > self.ymin) and (pt.y < self.ymax) and \
                    (pt.z > self.zmin) and (pt.z < self.zmax):
                return True
            else:
                return False
        except AttributeError:
            raise TypeError('point parameter must be a 3-element '
                            'sequence, ndarray, or Point object instance')


class Cube(Cuboid):
    """Class defining a cube.

    Parameters
    ----------
    a : float, optional
        length of edge
    center : {tuple, :class:`Point`}
        Either a 3-tuple of floats or an instance of the
        :class:`Point` class specifying the :math:`(x,y,z)` coordinates
        of the cube center.

    """
    def __init__(self, a=None, center=None):

        if a is None:
            raise TypeError('Please specify the edge length parameter `a`')

        self._a = a

        self._center = None
        if center is None:
            self._center = Point()
        elif isinstance(center, (tuple, list)):
            self._center = Point(x=center[0], y=center[1], z=center[2])
        elif isinstance(center, Point):
            self._center = center

        c0 = self._center
        cuboid_kwargs = {'xmin': c0.x - a / 2,
                         'ymin': c0.y - a / 2,
                         'zmin': c0.z - a / 2,
                         'xmax': c0.x + a / 2,
                         'ymax': c0.y + a / 2,
                         'zmax': c0.z + a / 2}

        super(Cube, self).__init__(**cuboid_kwargs)

    def __repr__(self):
        return("Cube(a={!r}, center={!r})".format(self.a, self.center))

    @property
    def a(self):
        return self._a

    @property
    def center(self):
        return self._center
