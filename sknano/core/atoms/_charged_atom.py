# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with a charge attribute (:mod:`sknano.core.atoms._charged_atom`)
===============================================================================

An `Atom` class with a charge attribute

.. currentmodule:: sknano.core.atoms._charged_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
import numbers
import numpy as np

from ._atom import Atom

__all__ = ['ChargedAtom']


@total_ordering
class ChargedAtom(Atom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    q : {int, float}, optional
        Net charge of `ChargedAtom`.
    """
    def __init__(self, *args, q=0., **kwargs):

        super().__init__(*args, **kwargs)

        self.q = q
        self.fmtstr = super().fmtstr + ", q={q!r}"

    def __eq__(self, other):
        return np.allclose(self.q, other.q) and super().__eq__(other)

    def __lt__(self, other):
        return (self.q < other.q and super().__le__(other)) or \
            (self.q <= other.q and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('q')
        return attrs

    @property
    def q(self):
        """Charge :math:`q` as multiple of elementary charge :math:`e`.

        """
        return self._q

    @q.setter
    def q(self, value):
        """Set `ChargedAtom` charge :math:`q`.

        Parameters
        ----------
        value : {int, float}
            net charge on atom as a multiple of the elementary charge
            :math:`e`.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._q = value

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(q=self.q))
        return super_dict
