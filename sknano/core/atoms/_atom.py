# -*- coding: utf-8 -*-
"""
============================================================================
Base class for structure data atom (:mod:`sknano.core.atoms._atom`)
============================================================================

.. currentmodule:: sknano.core.atoms._atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from functools import total_ordering

import numbers
import numpy as np

from sknano.core import xyz
from sknano.core.math import Vector, rotation_transform
from sknano.core.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_symbols, element_names

__all__ = ['Atom']


@total_ordering
class Atom(object):
    """Base class for structure data atom.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number :math:`\\boldsymbol{Z}`.
    x, y, z : float, optional
        :math:`x, y, z` coordinates of `Atom`.

    """
    _atomattrs = ['symbol', 'Z', 'm', 'r', 'x', 'y', 'z']

    def __init__(self, element=None, m=None, x=None, y=None, z=None):

        if m is None:
            m = 0

        if isinstance(element, numbers.Number):
            try:
                Z = int(element)
                idx = Z - 1
                symbol = element_symbols[idx]
                m = atomic_masses[symbol]
            except KeyError:
                print('unrecognized element number: {}'.format(element))
        else:
            if element in element_symbols:
                symbol = element
            elif element in element_names:
                symbol = element_symbols[element_names.index(element)]
            elif m in atomic_mass_symbol_map:
                symbol = atomic_mass_symbol_map[m]
            elif int(m / 2) in atomic_number_symbol_map:
                symbol = atomic_number_symbol_map[int(m / 2)]
            else:
                symbol = 'X'

        try:
            Z = atomic_numbers[symbol]
            m = atomic_masses[symbol]
        except KeyError:
            Z = 0

        self._symbol = symbol
        self._Z = Z
        self._m = m
        self._r = Vector([x, y, z])
        #self._p0 = Point([x, y, z])
        self._dr = Vector(np.zeros(3), p0=[x, y, z])

    def __str__(self):
        """Return string representation of atom."""
        return repr(self)

    def __repr__(self):
        """Return string representation of atom."""
        return "Atom(element={!r}, m={!r}, x={!r}, y={!r}, z={!r})".format(
            self.element, self.m, self.x, self.y, self.z)

    def __eq__(self, other):
        if self is other:
            return True
        else:
            for p in self._atomattrs:
                if getattr(self, p) != getattr(other, p):
                    return False
            return True

    def __lt__(self, other):
        return self.Z < other.Z

    @property
    def Z(self):
        """Atomic number :math:`Z`.

        Returns
        -------
        int
            Atomic number :math:`Z`.
        """
        return self._Z

    @property
    def element(self):
        """Element symbol.

        Returns
        -------
        str
            Symbol of chemical element.
        """
        return self.symbol

    @property
    def symbol(self):
        """Element symbol.

        Returns
        -------
        str
            Element symbol.
        """
        return self._symbol

    @property
    def m(self):
        """Atomic mass :math:`m_a` in atomic mass units.

        Returns
        -------
        float
            Atomic mass :math:`m_a` in atomic mass units.
        """
        return self._m

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
        self._r = Vector(value)

    @property
    def dr(self):
        """:math:`x, y, z` components of `Atom` displacement vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self._dr

    @dr.setter
    def dr(self, value):
        """Set :math:`x, y, z` components of `Atom` displacement vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` components of `Atom` displacement vector.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._dr = Vector(value)

    def todict(self):
        return dict(element=self.element, m=self.m,
                    x=self.x, y=self.y, z=self.z)

    def get_coords(self, asdict=False):
        """Return atom coords.

        Parameters
        ----------
        asdict : bool, optional

        Returns
        -------
        coords : :class:`python:~collections.OrderedDict` or ndarray

        """
        if asdict:
            return OrderedDict(zip(xyz, self.r))
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
        self._r.rezero(epsilon=epsilon)

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Atom.rezero`."""
        self.rezero(epsilon=epsilon)

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               deg2rad=False, transform_matrix=None):
        try:
            self.r = rotation_transform(self.r,
                                        transform_matrix=transform_matrix)
        except ValueError:
            self.r.rotate(angle, rot_axis=rot_axis, anchor_point=anchor_point,
                          deg2rad=deg2rad)

    def translate(self, t):
        #TODO compare timing benchmarks for translation of position vector.
        #self.r.translate(t, fix_tail=True)
        self.r += t
