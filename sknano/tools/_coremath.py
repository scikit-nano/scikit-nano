# -*- coding: utf-8 -*-
"""
==================================================================
Abstract data structures for math (:mod:`sknano.tools._coremath`)
==================================================================

.. currentmodule:: sknano.tools._coremath

"""
from __future__ import division, print_function, absolute_import

import numpy as np

from ._corefuncs import check_type

__all__ = ['Point', 'Vector', 'Quaternion']


class Point(object):
    """Create a point in :math:`R^3`

    Parameters
    ----------
    x, y, z : float, optional
        :math:`x, y, z` coordinates of point in :math:`R^3` space.
    units : {None, str}, optional
        Units of coordinates.

    """
    def __init__(self, x=None, y=None, z=None, units=None, dtype=None):
        self._p = np.zeros(3, dtype=float)
        self._units = units
        for i, pi in enumerate((x, y, z)):
            if pi is not None:
                self._p[i] = pi

    def __str__(self):
        return '({}, {}, {})'.format(self.x, self.y, self.z)

    def __repr__(self):
        return '({}, {}, {})'.format(self.x, self.y, self.z)

    @property
    def x(self):
        """:math:`x`-coordinate of `Point`.

        Returns
        -------
        float
            :math:`x`-coordinate of `Point`.

        """
        return self._p[0]

    @x.setter
    def x(self, value=float):
        """Set :math:`x`-coordinate of `Point`.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate of `Point`.

        """
        try:
            check_type(value, allowed_types=(int, float))
            self._p[0] = float(value)
        except TypeError as e:
            print(e)

    @property
    def y(self):
        """:math:`y`-coordinate of `Point`.

        Returns
        -------
        float
            :math:`y`-coordinate of `Point`.

        """
        return self._p[1]

    @y.setter
    def y(self, value=float):
        """Set :math:`y`-coordinate of `Point`.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate of `Point`.

        """
        try:
            check_type(value, allowed_types=(int, float))
            self._p[1] = float(value)
        except TypeError as e:
            print(e)

    @property
    def z(self):
        """:math:`z`-coordinate of `Point`.

        Returns
        -------
        float
            :math:`z`-coordinate of `Point`.

        """
        return self._p[2]

    @z.setter
    def z(self, value=float):
        """Set :math:`z`-coordinate of `Point`.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate of `Point`.

        """
        try:
            check_type(value, allowed_types=(int, float))
            self._p[2] = float(value)
        except TypeError as e:
            print(e)

    @property
    def coords(self):
        """:math:`x, y, z` coordinates of `Point`.

        Returns
        -------
        :py:class:`~numpy:numpy.ndarray`
            3-element :py:class:`~numpy:numpy.ndarray` of [:math:`x, y, z`]
            coordinates of `Point`.

        """
        return self._p

    @coords.setter
    def coords(self, value=np.ndarray):
        """Set :math:`x, y, z` coordinates of `Point`

        Parameters
        ----------
        value : :py:class:`~numpy:numpy.ndarray`
            3-element :py:class:`~numpy:numpy.ndarray` of [:math:`x, y, z`]
            coordinates of `Point`.

        """
        try:
            check_type(value, allowed_types=(np.ndarray,))
            for i, pi in enumerate(value):
                try:
                    check_type(pi, allowed_types=(int, float))
                    self._p[i] = pi
                except TypeError as e:
                    print(e)
        except TypeError as e:
            print(e)

    def fix_minus_zero_coords(self, epsilon=1.0e-10):
        """Set `Point` coordinates that are small, negative numbers to zero.

        Set `Point` coordinates that are negative and have absolute value
        less than `epsilon` to zero.

        Parameters
        ----------
        epsilon : float, optional
            Smallest allowed absolute value of any :math:`x,y,z` coordinate.

        """
        p = self._p.tolist()
        for i, pi in enumerate(p[:]):
            if pi < 0 and abs(pi) < epsilon:
                p[i] = 0.0
        self._p[0], self._p[1], self._p[2] = p

    def rezero_coords(self, epsilon=1.0e-10):
        """Re-zero `Point` coordinates near zero.

        Set `Point` coordinates with absolute value less than `epsilon` to
        zero.

        Parameters
        ----------
        epsilon : float, optional
            Smallest allowed absolute value of any :math:`x,y,z` coordinate.

        """
        p = self._p.tolist()
        for i, pi in enumerate(p[:]):
            if abs(pi) < epsilon:
                p[i] = 0.0
        self._p[0], self._p[1], self._p[2] = p


class Vector(object):
    """Create a vector in :math:`R^3`

    Parameters
    ----------
    x, y, z : float, optional
        :math:`x, y, z` components of terminating point of vector in
        :math:`R^3` space relative to origin.
    x0, y0, z0 : float, optional
        :math:`x_0, y_0, z_0` components of starting point of vector in
        :math:`R^3` space relative to origin.
    p, p0 : `Point`, optional
        Terminating and starting `Point` of vector in :math:`R^3` space
        relative to origin. If `p` is not `None` it will always override the
        `Point` defined by the `x`, `y`, `z` parameters. Similarly, if `p0` is
        not `None`, it will always override the `Point` defined by the
        `x0`, `y0`, `z0` parameters.
    units : {None, str}, optional
        Units of vector.


    Notes
    -----
    .. todo::

       add new methods for coordinate transformations


    """
    def __init__(self, x=None, y=None, z=None, x0=None, y0=None, z0=None,
                 p=None, p0=None, units=None):

        self._p = None
        if p is None:
            self._p = Point(x=x, y=y, z=z)
        elif isinstance(p, Point):
            self._p = p

        self._p0 = None
        if p0 is None:
            self._p0 = Point(x=x0, y=y0, z=z0)
        elif isinstance(p0, Point):
            self._p0 = p0

        self._v = np.zeros(3, dtype=float)
        for i, (pi, pi0) in enumerate(zip(self._p.coords, self._p0.coords)):
            self._v[i] = pi - pi0

    def __str__(self):
        return '({}, {}, {})'.format(self.x, self.y, self.z)

    def __repr__(self):
        return '({}, {}, {})'.format(self.x, self.y, self.z)

    @property
    def x(self):
        """:math:`x`-coordinate of `Vector`.

        Returns
        -------
        float
            :math:`x`-coordinate of `Vector`.

        """
        return self._v[0]

    @x.setter
    def x(self, value=float):
        """Set :math:`x`-coordinate of `Vector`.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate of `Vector`.

        """
        try:
            check_type(value, allowed_types=(int, float))
            self._v[0] = float(value)
        except TypeError as e:
            print(e)

    @property
    def y(self):
        """:math:`y`-coordinate of `Vector`.

        Returns
        -------
        float
            :math:`y`-coordinate of `Vector`.

        """
        return self._v[1]

    @y.setter
    def y(self, value=float):
        """Set :math:`y`-coordinate of `Vector`.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate of `Vector`.

        """
        try:
            check_type(value, allowed_types=(int, float))
            self._v[1] = float(value)
        except TypeError as e:
            print(e)

    @property
    def z(self):
        """:math:`z`-coordinate of `Vector`.

        Returns
        -------
        float
            :math:`z`-coordinate of `Vector`.

        """
        return self._v[2]

    @z.setter
    def z(self, value=float):
        """Set :math:`z`-coordinate of `Vector`.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate of `Vector`.

        """
        try:
            check_type(value, allowed_types=(int, float))
            self._v[2] = float(value)
        except TypeError as e:
            print(e)

    @property
    def components(self):
        """:math:`x, y, z` components of `Vector`.

        Returns
        -------
        :py:class:`~numpy:numpy.ndarray`
            3-element :py:class:`~numpy:numpy.ndarray` of [:math:`x, y, z`]
            components of `Vector`.

        """
        return self._v

    @components.setter
    def components(self, value=np.ndarray):
        """Set :math:`x, y, z` components of `Vector`

        Parameters
        ----------
        value : :py:class:`~numpy:numpy.ndarray`
            3-element :py:class:`~numpy:numpy.ndarray` of [:math:`x, y, z`]
            components of `Vector`.

        """
        try:
            check_type(value, allowed_types=(np.ndarray,))
            for i, vi in enumerate(value):
                try:
                    check_type(vi, allowed_types=(int, float))
                    self._v[i] = vi
                except TypeError as e:
                    print(e)
        except TypeError as e:
            print(e)

    def fix_minus_zero_components(self, epsilon=1.0e-10):
        """Set `Vector` components that are small, negative numbers to zero.

        Set `Vector` components that are negative and have absolute value
        less than `epsilon` to zero.

        Parameters
        ----------
        epsilon : float, optional
            Smallest allowed absolute value of any :math:`x,y,z` component.

        """
        v = self._v.tolist()
        for i, vi in enumerate(v[:]):
            if vi < 0 and abs(vi) < epsilon:
                v[i] = 0.0
        self._v[0], self._v[1], self._v[2] = v

    def rezero_components(self, epsilon=1.0e-10):
        """Re-zero `Vector` components near zero.

        Set `Vector` components with absolute value less than `epsilon` to
        zero.

        Parameters
        ----------
        epsilon : float, optional
            Smallest allowed absolute value of any :math:`x,y,z` component.

        """
        v = self._v.tolist()
        for i, vi in enumerate(v[:]):
            if abs(vi) < epsilon:
                v[i] = 0.0
        self._v[0], self._v[1], self._v[2] = v


class Quaternion(object):
    pass
