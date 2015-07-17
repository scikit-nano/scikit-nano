# -*- coding: utf-8 -*-
"""
============================================================================
Base class for atom bond (:mod:`sknano.core.atoms._bond`)
============================================================================

.. currentmodule:: sknano.core.atoms._bond

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering

import numpy as np

from sknano.core import BaseClass
from sknano.core.math import Vector

import sknano.core.atoms

__all__ = ['Bond']


@total_ordering
class Bond(BaseClass):
    """Abstract representation of bond between 2 `Atom` objects.

    Parameters
    ----------
    atom1, atom2 : `~sknano.core.atoms.Atom` instances

    """
    def __init__(self, atom1, atom2):
        super().__init__()
        self.atoms = sknano.core.atoms.StructureAtoms([atom1, atom2])
        self.vector = Vector(atom2.r - atom1.r, p0=atom1.r.p)
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
    def vector(self):
        """`Bond` :class:`~sknano.core.math.Vector`.

        `Bond` :class:`~sknano.core.math.Vector` points from
        :attr:`Bond.atom1` to :attr:`Bond.atom2`.

        """
        return self._vector

    @vector.setter
    def vector(self, value):
        """Set `Bond` :class:`~sknano.core.math.Vector`."""
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an `array_like` object')
        self._vector = Vector(value)

    @property
    def unit_vector(self):
        """`Bond` :attr:`~sknano.core.math.Vector.unit_vector`."""
        return self.vector.unit_vector

    @property
    def length(self):
        """`Bond` :attr:`~sknano.core.math.Vector.length`."""
        return self.vector.length

    def todict(self):
        return dict(atom1=self.atom1, atom2=self.atom2)
