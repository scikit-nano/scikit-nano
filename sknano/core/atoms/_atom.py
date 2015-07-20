# -*- coding: utf-8 -*-
"""
============================================================================
Base class for structure data atom (:mod:`sknano.core.atoms._atom`)
============================================================================

.. currentmodule:: sknano.core.atoms._atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering

import numbers

import numpy as np

from sknano.core import BaseClass
from sknano.core.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_symbols, element_names

__all__ = ['Atom']


@total_ordering
class Atom(BaseClass):
    """Base class for abstract representation of structure atom.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number :math:`\\boldsymbol{Z}`.

    """
    # _fields = ['element']

    def __init__(self, *args, element=None, mass=None, Z=None, **kwargs):
        args = list(args)

        if 'm' in kwargs and mass is None:
            mass = kwargs['m']
            del kwargs['m']

        if element is None:

            if mass is not None or Z is not None:
                if Z is not None:
                    args.append(Z)
                else:
                    args.append(mass)

            if len(args) > 0:
                element = args.pop()

        super().__init__(*args, **kwargs)
        self.mass = mass
        self.element = element
        self.fmtstr = "{element!r}, Z={Z!r}, mass={mass!r}"

    def __eq__(self, other):
        """Test equality of two `Atom` object instances."""
        return self is other or (self.element == other.element and
                                 self.Z == other.Z and
                                 np.allclose(self.mass, other.mass))

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        if self.element == other.element == 'X' and \
                self.Z == other.Z == 0:
            return self.mass < other.mass
        else:
            return self.Z < other.Z

    def __dir__(self):
        return ['element', 'Z', 'mass']

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
        symbol = None
        if isinstance(value, numbers.Integral):
            try:
                Z = int(value)
                idx = Z - 1
                symbol = element_symbols[idx]
            except IndexError:
                print('unrecognized element number: {}'.format(value))

        if symbol is None:
            if isinstance(value, str):
                if value in element_symbols:
                    symbol = value
                elif value.capitalize() in element_names:
                    symbol = element_symbols[element_names.index(value)]
            elif isinstance(value, numbers.Number):
                if value in atomic_mass_symbol_map:
                    symbol = atomic_mass_symbol_map[value]
                elif int(value / 2) in atomic_number_symbol_map:
                    symbol = atomic_number_symbol_map[int(value / 2)]

        if symbol is None:
            symbol = 'X'

        try:
            self._Z = atomic_numbers[symbol]
            self._mass = atomic_masses[symbol]
        except KeyError:
            self._Z = 0
            if self.mass is None:
                self._mass = 0

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

    def rezero(self, *args, **kwargs):
        assert not hasattr(super(), 'rezero')

    def rotate(self, **kwargs):
        assert not hasattr(super(), 'rotate')

    def translate(self, *args, **kwargs):
        assert not hasattr(super(), 'translate')

    def todict(self):
        """Return `dict` of `Atom` constructor parameters."""
        return dict(element=self.element, mass=self.mass, Z=self.Z)
