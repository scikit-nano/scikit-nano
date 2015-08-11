# -*- coding: utf-8 -*-
"""
===============================================================================
Class representation of atom bonds (:mod:`sknano.core.atoms._bonds`)
===============================================================================

.. currentmodule:: sknano.core.atoms._bonds

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter

import numpy as np

from sknano.core import BaseClass, UserList, cyclic_pairs, dedupe
from sknano.core.math import Vector, vector as vec

import sknano.core.atoms

__all__ = ['Bond', 'Bonds']


@total_ordering
class Bond(BaseClass):
    """Abstract representation of bond between 2 `Atom` objects.

    Parameters
    ----------
    atom1, atom2 : `~sknano.core.atoms.Atom` instances

    """
    def __init__(self, atom1, atom2):
        super().__init__()
        self.atoms = sknano.core.atoms.StructureAtoms()
        self.atoms.extend([atom1, atom2])
        self.fmtstr = "{atom1!r}, {atom2!r}"

    def __str__(self):
        """Return nice string representation of `Bond`."""
        return "Bond({!r}->{!r})".format(self.atom1.id, self.atom2.id)

    def __eq__(self, other):
        return self is other or (self.atoms == other.atoms)

    def __lt__(self, other):
        return self.length < other.length

    def __dir__(self):
        return ['atoms', 'atom1', 'atom2', 'vector', 'unit_vector', 'length']

    @property
    def atoms(self):
        """:class:`~sknano.core.atoms.Atoms` in `Bond`."""
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        self._atoms = value

    @property
    def atom1(self):
        """:class:`~sknano.core.atoms.Atom` 1 in `Bond`."""
        return self.atoms[0]

    @property
    def atom2(self):
        """:class:`~sknano.core.atoms.Atom` 2 in `Bond`."""
        return self.atoms[1]

    @property
    def centroid(self):
        """:attr:`~sknano.core.atoms.XYZAtoms.centroid` of :class:`Bond` \
            :attr:`~Bond.atoms`."""
        return self.atoms.centroid

    @property
    def vector(self):
        """`Bond` :class:`~sknano.core.math.Vector`.

        `Bond` :class:`~sknano.core.math.Vector` points from
        :attr:`Bond.atom1` to :attr:`Bond.atom2`.

        """
        return Vector(self.atom2.r - self.atom1.r, p0=self.atom1.r.p)

    @property
    def unit_vector(self):
        """`Bond` :attr:`~sknano.core.math.Vector.unit_vector`."""
        return self.vector.unit_vector

    @property
    def length(self):
        """`Bond` :attr:`~sknano.core.math.Vector.length`."""
        return self.vector.length

    def rotate(self, **kwargs):
        """Rotate the bond by rotating the :attr:`~Bond.atoms`."""
        [atom.rotate(fix_anchor_point=True, **kwargs) for atom in self.atoms]

    def todict(self):
        return dict(atom1=self.atom1, atom2=self.atom2)


class Bonds(UserList):
    """Base class for collection of atom `Bond`\ s.

    Parameters
    ----------
    bonds : {None, sequence, `Bonds`}, optional
        if not `None`, then a list of `Bond` instance objects

    """
    def __init__(self, bonds=None):
        super().__init__(initlist=bonds)
        self.fmtstr = "{bonds!r}"

    def __repr__(self):
        """Return canonical string representation of `Bonds`."""
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def sort(self, key=attrgetter('length'), reverse=False):
        super().sort(key=key, reverse=reverse)

    # def __getitem__(self, index):
    #     data = super().__getitem__(index)
    #     if isinstance(data, list):
    #         return self.__class__(data)
    #     return data

    # def __setitem__(self, index, value):
    #     if not isinstance(value, (self.__class__, Bond)):
    #         if isinstance(index, slice):
    #             value = self.__class__(value)
    #         else:
    #             value = Bond(value)
    #     super().__setitem__(index, value)

    @property
    def fmtstr(self):
        return self._fmtstr

    @fmtstr.setter
    def fmtstr(self, value):
        self._fmtstr = value

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
        return cyclic_pairs(self.data)

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
        atoms = dedupe(atoms, key=attrgetter('id'))
        return sknano.core.atoms.StructureAtoms(list(atoms))

    def todict(self):
        return dict(bonds=self.data)
