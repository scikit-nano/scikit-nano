# -*- coding: utf-8 -*-
"""
============================================================================
Base class for structure data atom (:mod:`sknano.core.atoms._atom`)
============================================================================

.. currentmodule:: sknano.core.atoms._atom

"""
from __future__ import absolute_import, division, print_function

__docformat__ = 'restructuredtext en'

from functools import total_ordering

import numbers

from sknano.core.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_symbols, element_names

__all__ = ['Atom']


@total_ordering
class Atom(object):
    """Base class for abstract representation of structure atom.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number :math:`\\boldsymbol{Z}`.

    """
    _atomattrs = ['symbol', 'Z', 'mass']  # private class var

    __hash__ = object.__hash__

    def __init__(self, element=None, mass=None, **kwargs):

        # set mass first because the element.setter method may check mass value
        if mass is None and 'm' in kwargs:
            mass = kwargs['m']

        if mass is None:
            mass = 0
        self.mass = mass

        self.element = element

    def __str__(self):
        """Return nice string representation of `Atom`."""
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `Atom`."""
        return "Atom(element={!r}, mass={!r})".format(self.element, self.mass)

    def __eq__(self, other):
        """Test equality of two `Atom` object instances."""
        if self is other:
            return True
        else:
            for p in self._atomattrs:
                if getattr(self, p) != getattr(other, p):
                    return False
            return True

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        return self.Z < other.Z

    #def __mul__(self, other):
    #    """Multiply atom

    @property
    def Z(self):
        """Atomic number :math:`Z`.

        Returns
        -------
        int
            Atomic number :math:`Z`.
        """
        return self._Z

    @Z.setter
    def Z(self, value):
        """Set atomic number :math:`Z`.

        Parameters
        ----------
        value : int
            Atomic number :math:`Z`.

        """
        if not (isinstance(value, numbers.Real) and int(value) > 0):
            raise ValueError('Expected a real, positive integer.')
        try:
            Z = int(value)
            idx = Z - 1
            symbol = element_symbols[idx]
            mass = atomic_masses[symbol]
        except KeyError:
            print('unrecognized element number: {}'.format(value))
        else:
            self._Z = atomic_numbers[symbol]
            self._mass = mass
            self._symbol = symbol

    @property
    def element(self):
        """Element symbol.

        Returns
        -------
        str
            Symbol of chemical element.
        """
        return self._symbol

    @element.setter
    def element(self, value):
        """Set element symbol."""
        if isinstance(value, numbers.Number):
            try:
                Z = int(value)
                idx = Z - 1
                symbol = element_symbols[idx]
                mass = atomic_masses[symbol]
            except KeyError:
                print('unrecognized element number: {}'.format(value))
            else:
                self._Z = atomic_numbers[symbol]
                self._mass = mass
        else:
            mass = self.mass
            if value in element_symbols:
                symbol = value
            elif value in element_names:
                symbol = element_symbols[element_names.index(value)]
            elif mass in atomic_mass_symbol_map:
                symbol = atomic_mass_symbol_map[mass]
            elif int(mass / 2) in atomic_number_symbol_map:
                symbol = atomic_number_symbol_map[int(mass / 2)]
            else:
                symbol = 'X'

            try:
                self._Z = atomic_numbers[symbol]
                self._mass = atomic_masses[symbol]
            except KeyError:
                self._Z = 0

        self._symbol = symbol

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
    def mass(self):
        """Atomic mass :math:`m_a` in atomic mass units.

        Returns
        -------
        float
            Atomic mass :math:`m_a` in atomic mass units.
        """
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value

    @property
    def m(self):
        return self.mass

    @m.setter
    def m(self, value):
        self.mass = value

    def todict(self):
        """Return `dict` of `Atom` constructor parameters."""
        return dict(element=self.element, mass=self.mass)
