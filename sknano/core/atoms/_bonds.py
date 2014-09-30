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

import numpy as np

from sknano.core.math import vector as vec
from ._extended_atoms import XAtoms
from ._base import UserList
#from ._bond import Bond

__all__ = ['Bonds']


class Bonds(UserList):
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
        super(Bonds, self).__init__(initlist=bonds, copylist=copylist,
                                    deepcopy=deepcopy)

    def __str__(self):
        """Return a nice string representation of `Bonds`."""
        return "Bonds({!s})".format(self.data)

    def __repr__(self):
        """Return the canonical string representation of `Bonds`."""
        return "Bonds({!r})".format(self.data)

    def sort(self, key=None, reverse=False):

        if key is None:
            self.data.sort(key=attrgetter('length'), reverse=reverse)
        else:
            self.data.sort(key=key, reverse=reverse)

    @property
    def Nbonds(self):
        """Number of `Bond`s in `Bonds`."""
        return len(self)

    @property
    def vectors(self):
        return np.asarray([bond.vector for bond in self])

    @property
    def lengths(self):
        return np.asarray([bond.length for bond in self])

    @property
    def mean_length(self):
        return np.mean(self.lengths)

    @property
    def angles(self):
        return np.asarray([vec.angle(b1.vector, b2.vector) for (b1, b2) in
                           combinations(self, 2)])

    @property
    def mean_angle(self):
        return np.mean(self.angles)

    @property
    def atoms(self):
        atoms = XAtoms()
        [atoms.extend(bond.atoms) for bond in self]
        return XAtoms(atoms=list(set(atoms)))
