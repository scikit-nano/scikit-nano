# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with energy attributes (:mod:`sknano.core.atoms._energy_atom`)
===============================================================================

An `Atom` class with energy attributes

.. currentmodule:: sknano.core.atoms._energy_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import zip

__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from functools import total_ordering
import numbers
import numpy as np

from sknano.core import xyz
from sknano.core.math import Point, Vector
from ._atom import Atom

__all__ = ['EnergyAtom']


@total_ordering
class EnergyAtom(Atom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.

    """
    def __init__(self, *args, pe=None, ke=None, etotal=None, **kwargs):

        super().__init__(*args, **kwargs)

        self.pe = pe
        self.ke = ke
        self.etotal = etotal

        self.fmtstr = super().fmtstr + \
            ", pe={pe!r}, ke={ke!r}, etotal={etotal!r}"

    def __eq__(self, other):
        return super().__eq__(other)

    def __lt__(self, other):
        if self.etotal < other.etotal:
            return True
        else:
            return super().__lt__(other)

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['pe', 'ke', 'etotal'])
        return attrs

    @property
    def pe(self):
        """`Atom` potential energy."""
        return self._pe

    @pe.setter
    def pe(self, value):
        if not isinstance(value, (float, type(None))):
            raise TypeError('Expected a number')
        if value is None:
            value = 0.0
        self._pe = value

    @property
    def ke(self):
        """`Atom` kinetic energy."""
        return self._ke

    @ke.setter
    def ke(self, value):
        if not isinstance(value, (float, type(None))):
            raise TypeError('Expected a number')
        if value is None:
            value = 0.0
        self._ke = value

    @property
    def etotal(self):
        """`Atom` total energy \
                (:attr:`~EnergyAtom.pe` + :attr:`~EnergyAtom.ke`)."""
        return self._etotal

    @etotal.setter
    def etotal(self, value):
        if not isinstance(value, (float, type(None))):
            raise TypeError('Expected a number')
        if value is None:
            try:
                value = self.pe + self.ke
            except TypeError:
                value = 0.0
        self._etotal = value

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(pe=self.pe, ke=self.ke, etotal=self.etotal))
        return super_dict
