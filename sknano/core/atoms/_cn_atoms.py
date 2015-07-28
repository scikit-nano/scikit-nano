# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for `CNAtom`\ s (:mod:`sknano.core.atoms._cn_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._cn_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Counter
from operator import attrgetter

import numpy as np

from ._atoms import Atoms
from ._cn_atom import CNAtom

__all__ = ['CNAtoms']


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
