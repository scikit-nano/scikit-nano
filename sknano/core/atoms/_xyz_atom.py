# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with x, y, z attributes (:mod:`sknano.core.atoms._xyz_atom`)
===============================================================================

An `Atom` sub-class with x, y, z attributes

.. currentmodule:: sknano.core.atoms._xyz_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from functools import total_ordering
import numbers
import numpy as np

from sknano.core import xyz
from sknano.core.math import Vector
from ._atom import Atom

__all__ = ['XYZAtom']


@total_ordering
class XYZAtom(Atom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    x, y, z : float, optional
        :math:`x, y, z` components of `XYZAtom` position vector relative to
        origin.

    """
    def __init__(self, *args, x=None, y=None, z=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.r0 = Vector([x, y, z])
        self.r = Vector([x, y, z])
        self.fmtstr = super().fmtstr + ", x={x:.6f}, y={y:.6f}, z={z:.6f}"

    def __eq__(self, other):
        return self.r == other.r and super().__eq__(other)

    def __lt__(self, other):
        return (self.r < other.r and super().__le__(other)) or \
            (self.r <= other.r and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['x', 'y', 'z'])
        return attrs

    @property
    def x(self):
        """:math:`x`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        return self._r.x

    @x.setter
    def x(self, value):
        """Set `Atom` :math:`x`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.x = value

    @property
    def y(self):
        """:math:`y`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        return self._r.y

    @y.setter
    def y(self, value):
        """Set `Atom` :math:`y`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.y = value

    @property
    def z(self):
        """:math:`z`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        return self._r.z

    @z.setter
    def z(self, value):
        """Set `Atom` :math:`z`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.z = value

    @property
    def r(self):
        """:math:`x, y, z` components of `Atom` position vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self._r

    @r.setter
    def r(self, value):
        """Set :math:`x, y, z` components of `Atom` position vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Atom` position vector relative to
            the origin.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._r = Vector(value, nd=3)

    @property
    def r0(self):
        """:math:`x, y, z` components of `Atom` position vector at t=0.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self._r0

    @r0.setter
    def r0(self, value):
        """Set :math:`x, y, z` components of `Atom` position vector at t=0.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Atom` position vector relative to
            the origin.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._r0 = Vector(value, nd=3)

    @property
    def dr(self):
        """:math:`x, y, z` components of `Atom` displacement vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self.r - self.r0

    def get_coords(self, asdict=False):
        """Return atom coords.

        Parameters
        ----------
        asdict : bool, optional

        Returns
        -------
        coords : :class:`~python:collections.OrderedDict` or ndarray

        """
        if asdict:
            return OrderedDict(list(zip(xyz, self.r)))
        else:
            return self.r

    def rezero(self, epsilon=1.0e-10):
        """Re-zero position vector components.

        Set position vector components with absolute value less than
        `epsilon` to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        self._r.rezero(epsilon)
        super().rezero(epsilon)

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Atom.rezero`."""
        self.rezero(epsilon=epsilon)

    def rotate(self, **kwargs):
        """Rotate `Atom` position vector.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        self.r.rotate(**kwargs)
        super().rotate(**kwargs)

    def translate(self, t, fix_anchor_point=True):
        """Translate `Atom` position vector by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : bool, optional

        """
        # TODO compare timing benchmarks for translation of position vector.
        self.r.translate(t, fix_anchor_point=fix_anchor_point)
        super().translate(t, fix_anchor_point=fix_anchor_point)
        # self.r += t

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(x=self.x, y=self.y, z=self.z))
        return super_dict
