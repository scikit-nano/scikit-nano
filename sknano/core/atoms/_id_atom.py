# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with extended feature set (:mod:`sknano.core.atoms._extended_atom`)
===============================================================================

An "eXtended" `Atom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._extended_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
import numbers
# import numpy as np

from ._atom import Atom

__all__ = ['IDAtom']


@total_ordering
class IDAtom(Atom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    id : int, optional
        atom ID
    mol : int, optional
        molecule ID
    """
    def __init__(self, *args, id=0, mol=0, **kwargs):

        if 'atomID' in kwargs:
            id = kwargs['atomID']
            del kwargs['atomID']

        if 'moleculeID' in kwargs:
            mol = kwargs['moleculeID']
            del kwargs['moleculeID']

        super().__init__(*args, **kwargs)

        self.id = id
        self.mol = mol
        self.fmtstr = super().fmtstr + ", id={id!r}, mol={mol!r}"

    def __eq__(self, other):
        return self.id == other.id and self.mol == other.mol and \
            super().__eq__(other)

    def __lt__(self, other):
        return ((self.id < other.id and self.mol <= other.mol and
                 super().__le__(other))
                or (self.id <= other.id and self.mol < other.mol and
                    super().__le__(other))
                or (self.id <= other.id and self.mol <= other.mol and
                    super().__lt__(other)))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['id', 'mol'])
        return attrs

    @property
    def id(self):
        """:attr:`~IDAtom.id`."""
        return self._id

    @id.setter
    def id(self, value):
        """Set atom id.

        Parameters
        ----------
        value : int
            atom ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._id = int(value)

    @property
    def mol(self):
        """:attr:`~IDAtom.mol`."""
        return self._mol

    @mol.setter
    def mol(self, value):
        """Set :attr:`~IDAtom.mol`.

        Parameters
        ----------
        value : int
            molecule ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._mol = int(value)

    @property
    def atomID(self):
        return self.id

    @atomID.setter
    def atomID(self, value):
        self.id = value

    @property
    def moleculeID(self):
        return self.mol

    @moleculeID.setter
    def moleculeID(self, value):
        self.mol = value

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(id=self.id, mol=self.mol))
        return super_dict
