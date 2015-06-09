# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with a coordination number (:mod:`sknano.core.atoms._cn_atom`)
===============================================================================

.. currentmodule:: sknano.core.atoms._cn_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
import numbers

import numpy as np

from ._atom import Atom

__all__ = ['CNAtom']


@total_ordering
class CNAtom(Atom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    CN : {int}, optional
        Coordination number.

    """
    def __init__(self, *args, CN=0, **kwargs):

        super().__init__(*args, **kwargs)

        self.CN = CN
        self.fmtstr = super().fmtstr + ", CN={CN!r}"

    def __eq__(self, other):
        return np.allclose(self.CN, other.CN) and super().__eq__(other)

    def __lt__(self, other):
        return (self.CN < other.CN and super().__le__(other)) or \
            (self.CN <= other.CN and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('CN')
        return attrs

    @property
    def CN(self):
        """Return `CNAtom` coordination number."""
        return self._CN

    @CN.setter
    def CN(self, value):
        """Set `CNAtom` coordination number."""
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number.')
        self._CN = int(value)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(CN=self.CN))
        return super_dict
