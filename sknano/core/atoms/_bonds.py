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

from operator import attrgetter

import numpy as np

from sknano.core import BaseClass, UserList, cyclic_pairs, dedupe, flatten, \
    pairwise
from sknano.core.math import Vector, vector as vec

import sknano.core.atoms

__all__ = ['Bond', 'Bonds']


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

    def _is_valid_operand(self, other):
        return isinstance(other, self.__class__)

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.atoms == other.atoms and super().__eq__(other)

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if not super().__le__(other) or self.atoms > other.atoms:
            return False
        return True

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return (self.atoms < other.atoms and self.__le__(other))

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if not super().__ge__(other) or self.atoms < other.atoms:
            return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.atoms > other.atoms and self.__ge__(other)

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

    @property
    def __item_class__(self):
        return Bond

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
        return list(cyclic_pairs(self.data)) if self.Nbonds > 2 \
            else (list(pairwise(self.data)) if self.Nbonds > 1 else None)

    @property
    def angles(self):
        """:class:`~numpy:numpy.ndarray` of `Bond` pair angles."""
        if self.Nbonds > 2:
            try:
                return np.asarray([vec.angle(b1.vector, b2.vector) for
                                   (b1, b2) in cyclic_pairs(self)])
            except TypeError:
                return None
        elif self.Nbonds == 2:
            return vec.angle(self[0].vector, self[-1].vector)
        else:
            return None

    @property
    def mean_angle(self):
        """Mean bond angle."""
        try:
            return np.mean(self.angles)
        except TypeError:
            return None

    @property
    def atoms(self):
        """`Atoms` :class:`python:set` in `Bonds`."""
        atoms = []
        [atoms.extend(bond.atoms) for bond in self]
        atoms = list(dedupe(list(flatten([bond.atoms for bond in self])),
                            key=attrgetter('id')))
        return sknano.core.atoms.StructureAtoms(atoms)

    def todict(self):
        return dict(bonds=self.data)
