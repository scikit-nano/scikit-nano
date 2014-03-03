# -*- coding: utf-8 -*-
"""
============================================================================
Base class for structure data atom (:mod:`sknano.structure_io.atoms._atom`)
============================================================================

.. currentmodule:: sknano.structure_io.atoms._atom

"""
from __future__ import division, absolute_import, print_function
__docformat__ = 'restructuredtext'

import numpy as np

from ...tools.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_symbols

__all__ = ['Atom', 'AtomAttributes']


class AtomAttributes(object):
    pass


class Atom(object):
    """Base class for structure data atom.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    x, y, z : float, optional
        :math:`x, y, z` components of `Atom` position vector relative to
        origin.

    """

    def __init__(self, element=None, mass=None, x=None, y=None, z=None):

        self._symbol = None
        self._Z = None
        self._m = None

        self._r = np.zeros(3, dtype=float)
        for i, ri in enumerate((x, y, z)):
            if ri is not None:
                self._r[i] = ri

        if isinstance(element, (int, float)):
            self._Z = int(element)
            idx = self._Z - 1
            try:
                self._symbol = element_symbols[idx]
                self._m = atomic_masses[self._symbol]
            except KeyError as e:
                print(e)
                print('unrecognized element number: {}'.format(element))
        elif isinstance(element, str):
            self._symbol = element
            try:
                self._Z = atomic_numbers[self._symbol]
                self._m = atomic_masses[self._symbol]
            except KeyError as e:
                print(e)
                print('Unrecognized atomic symbol: {}'.format(element))
        else:
            self._symbol = None
            self._Z = None
            if mass is not None and isinstance(mass, (int, float)):
                try:
                    if isinstance(mass, float):
                        self._symbol = atomic_mass_symbol_map[mass]
                    elif isinstance(mass, int):
                        self._symbol = atomic_number_symbol_map[int(mass / 2)]
                    self._Z = atomic_numbers[self._symbol]
                    self._m = atomic_masses[self._symbol]
                except KeyError as e:
                    self._symbol = None
                    self._Z = None
                    self._m = mass
            else:
                self._m = 0

        self._attributes = ['symbol', 'Z', 'm', 'r']

    def __str__(self):
        """Return string representation of atom."""
        atom_str = ''
        for attr in self._attributes:
            atom_str += \
                'Atom {}: {}\n'.format(attr, getattr(self, '_' + attr))
        return atom_str

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
        return self._r[0]

    @x.setter
    def x(self, value=float):
        """Set `Atom` :math:`x`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        self._check_type(value, (int, float))
        self._r[0] = float(value)

    @property
    def y(self):
        """:math:`y`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        return self._r[1]

    @y.setter
    def y(self, value=float):
        """Set `Atom` :math:`y`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        self._check_type(value, (int, float))
        self._r[1] = float(value)

    @property
    def z(self):
        """:math:`z`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        return self._r[2]

    @z.setter
    def z(self, value=float):
        """Set `Atom` :math:`z`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        self._check_type(value, (int, float))
        self._r[2] = float(value)

    @property
    def r(self):
        """:math:`x, y, z` position in units of **Angstroms**.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates
            of `Atom`.

        """
        return self._r

    @r.setter
    def r(self, value=np.ndarray):
        """Set :math:`x, y, z` coordinates of atom.

        Parameters
        ----------
        value : array_like
            3-element array of :math:`x, y, z`-coordinates in units of
            **Angstroms**.

        """
        self._check_type(value, np.ndarray)
        for i, ri in enumerate(value):
            self._r[i] = ri

    def rezero_coords(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        r = self._r.tolist()
        for i, ri in enumerate(r[:]):
            if abs(ri) < epsilon:
                r[i] = 0.0
        self._r[0], self._r[1], self._r[2] = r
