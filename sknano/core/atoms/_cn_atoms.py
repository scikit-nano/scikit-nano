# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with a coordination number (:mod:`sknano.core.atoms._cn_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._cn_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from collections import Counter
from functools import total_ordering
from operator import attrgetter
import numbers

import numpy as np

from ._atoms import Atom, Atoms

__all__ = ['CNAtom', 'CNAtoms']


@total_ordering
class CNAtom(Atom):
    """An `Atom` class with a coordination number attribute.

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


class CNAtoms(Atoms):
    """An `Atoms` sub-class for `CNAtom`\ s.

    A container class for :class:`~sknano.core.atoms.CNAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `CNAtoms`}, optional
        if not `None`, then a list of `CNAtom` instance objects or an
        existing `CNAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return CNAtom

    def sort(self, key=attrgetter('CN'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def coordination_numbers(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`CNAtom.CN`\ s."""
        return np.asarray([atom.CN for atom in self])

    @property
    def coordination_number_counts(self):
        """Coordination number counts.

        Returns
        -------
        :class:`~python:collections.Counter`
            :class:`~python:collections.Counter` of
            :attr:`~CNAtoms.coordination_numbers`.

        """
        return Counter(self.coordination_numbers)

    @property
    def coordination_counts(self):
        """Alias for :attr:`~CNAtoms.coordination_number_counts`."""
        return self.coordination_number_counts

    @property
    def CN_counts(self):
        """Alias for :attr:`~CNAtoms.coordination_number_counts`."""
        return self.coordination_number_counts
