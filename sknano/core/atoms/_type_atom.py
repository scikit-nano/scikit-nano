# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with a 'type' attribute (:mod:`sknano.core.atoms._type_atom`)
===============================================================================

An `Atom` class with a `type` attribute.

.. currentmodule:: sknano.core.atoms._type_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
import numbers

from ._atom import Atom

__all__ = ['TypeAtom']


@total_ordering
class TypeAtom(Atom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    type : int, optional
        atom type

    """
    def __init__(self, *args, type=1, **kwargs):
        if 'atomtype' in kwargs:
            type = kwargs['atomtype']
            del kwargs['atomtype']

        super().__init__(*args, **kwargs)

        self.type = type
        self.fmtstr = super().fmtstr + ", type={type!r}"

    def __eq__(self, other):
        return self.type == other.type and super().__eq__(other)

    def __lt__(self, other):
        return (self.type < other.type and super().__le__(other)) or \
            (self.type <= other.type and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('type')
        return attrs

    @property
    def type(self):
        """:attr:`~TypeAtom.type`."""
        return self._type

    @type.setter
    def type(self, value):
        """Set :attr:`~TypeAtom.type`.

        Parameters
        ----------
        value : int
            atom type

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._type = int(value)

    @property
    def atomtype(self):
        return self.type

    @atomtype.setter
    def atomtype(self, value):
        self.type = value

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(type=self.type))
        return super_dict
