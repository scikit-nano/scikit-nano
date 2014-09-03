# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for collection of `Bond`s (:mod:`sknano.core.atoms._bonds`)
===============================================================================

.. currentmodule:: sknano.core.atoms._bonds

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from itertools import combinations
from operator import attrgetter
import copy

import numpy as np

from sknano.core.math import vector as vec
from ._base import BondList
#from ._bond import Bond

__all__ = ['Bonds']


class Bonds(BondList):
    """Base class for collection of atom `Bonds`.

    Parameters
    ----------
    bonds : {None, sequence, `Bonds`}, optional
        if not `None`, then a list of `Bond` instance objects
    copylist : bool, optional
        perform shallow copy of bonds list
    deepcopy : bool, optional
        perform deepcopy of bonds list

    """
    def __init__(self, bonds=None, copylist=True, deepcopy=False):
        self._data = []
        self._angles = []

        if bonds is not None:
            try:
                if copylist and not deepcopy:
                    self._data.extend(bonds[:])
                elif deepcopy:
                    self._data.extend(copy.deepcopy(bonds))
                else:
                    self._data.extend(bonds)
            except AttributeError:
                raise TypeError('`bonds={!r}` '.format(bonds) +
                                'is not a valid `Bonds` constructor '
                                'argument.\n bonds must be `None`, a list '
                                'of `Bond` objects, or a `Bonds` object '
                                'instance.')

    def __str__(self):
        """Return a nice string representation of `Bonds`."""
        return "Bonds({!s})".format(self._data)

    def __repr__(self):
        """Return the canonical string representation of `Bonds`."""
        return "Bonds({!r})".format(self._data)

    def sort(self, key=None, reverse=False):

        if key is None:
            self._data.sort(key=attrgetter('length'), reverse=reverse)
        else:
            self._data.sort(key=key, reverse=reverse)

    @property
    def Nbonds(self):
        return len(self)

    @property
    def vectors(self):
        return np.asarray([bond.vector for bond in self])

    @property
    def lengths(self):
        return np.asarray([bond.length for bond in self])

    @property
    def angles(self):
        self._compute_bond_angles()
        return np.asarray(self._angles)

    def _compute_bond_angles(self):
        try:
            self._angles = [vec.angle(b1.vector, b2.vector) for (b1, b2) in
                            combinations(self, 2)]
        except ValueError:
            pass
