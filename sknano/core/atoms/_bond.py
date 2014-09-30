# -*- coding: utf-8 -*-
"""
============================================================================
Base class for atom bond (:mod:`sknano.core.atoms._bond`)
============================================================================

.. currentmodule:: sknano.core.atoms._bond

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from functools import total_ordering

#import numpy as np

from sknano.core.math import Vector
from ._extended_atoms import XAtoms

__all__ = ['Bond']


@total_ordering
class Bond(object):
    """Base class for `Atom`-`Atom` bond.

    Parameters
    ----------

    """

    def __init__(self, atom1, atom2):

        self._vector = Vector(atom2.r - atom1.r, p0=atom1.r.p)
        self._atoms = XAtoms(atoms=[atom1, atom2])

    def __str__(self):
        """Return a nice string representation of atom."""
        strrep = "{!s}={!s}\nv={!s}\n|v|={!s}"
        return strrep.format(self.atom1.element, self.atom2.element,
                             self.vector, self.length)

    def __repr__(self):
        """Return canonical string representation of atom."""
        return "Bond({!r}, {!r})".format(self.atom1, self.atom2)

    def __eq__(self, other):
        if self is other or self.length == other.length:
            return True
        else:
            return False

    def __lt__(self, other):
        return self.length < other.length

    @property
    def atoms(self):
        """:class:`~sknano.core.atoms.XAtoms` in `Bond`."""
        return self._atoms

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
        self._vector = Vector(value)

    @property
    def length(self):
        """`Bond` :attr:`~sknano.core.math.Vector.length`."""
        return self.vector.length
