# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for collection of `Bond`\ s (:mod:`sknano.core.atoms._bonds`)
===============================================================================

.. currentmodule:: sknano.core.atoms._bonds

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import str
__docformat__ = 'restructuredtext en'

#from itertools import combinations
from operator import attrgetter

import numpy as np

from sknano.core import UserList, cyclic_pairs, dedupe
from sknano.core.math import vector as vec

import sknano.core.atoms
#from ._bond import Bond

__all__ = ['Bonds']


class Bonds(UserList):
    """Base class for collection of atom `Bond`\ s.

    Parameters
    ----------
    bonds : {None, sequence, `Bonds`}, optional
        if not `None`, then a list of `Bond` instance objects

    """
    def __init__(self, bonds=None):
        super().__init__(initlist=bonds)

    def __str__(self):
        """Return a nice string representation of `Bonds`."""
        return '\n'.join([str(bond) for bond in self])

    def __repr__(self):
        """Return the canonical string representation of `Bonds`."""
        return "Bonds({!r})".format(self.data)

    def sort(self, key=attrgetter('length'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def Nbonds(self):
        """Number of `Bond`\ s in `Bonds`."""
        return len(self)

    @property
    def vectors(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Bond.vector`\ s."""
        return np.asarray([bond.vector for bond in self])

    @property
    def unit_vectors(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Bond.unit_vector`\ s."""
        return np.asarray([bond.unit_vector for bond in self])

    @property
    def lengths(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Bond.length`\ s."""
        return np.asarray([bond.length for bond in self])

    @property
    def mean_length(self):
        """Mean bond length."""
        return np.mean(self.lengths)

    @property
    def bond_angle_pairs(self):
        """`cyclic_pairs` of `Bond`\ s."""
        return cyclic_pairs(self)

    @property
    def angles(self):
        """:class:`~numpy:numpy.ndarray` of `Bond` pair angles."""
        return np.asarray([vec.angle(b1.vector, b2.vector) for (b1, b2) in
                           cyclic_pairs(self)])

    @property
    def mean_angle(self):
        """Mean bond angle."""
        return np.mean(self.angles)

    @property
    def atoms(self):
        """`Atoms` :class:`python:set` in `Bonds`."""
        atoms = []
        [atoms.extend(bond.atoms) for bond in self]
        atoms = dedupe(atoms, key=attrgetter('id', 'x', 'y', 'z'))
        return sknano.core.atoms.StructureAtoms(atoms=list(atoms))
