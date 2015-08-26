# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with bonds (:mod:`sknano.core.atoms._bond_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._bond_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter

import numpy as np

from ._atoms import Atom, Atoms
from ._bonds import Bond, Bonds

__all__ = ['BondAtom', 'BondAtoms']


@total_ordering
class BondAtom(Atom):
    """An `Atom` class for atom bonds.

    Parameters
    ----------
    bonds : {int}, optional
        Coordination number.

    """
    def __init__(self, *args, bonded_atoms=None, **kwargs):

        super().__init__(*args, **kwargs)

        self.bonds = [()]
        self.fmtstr = super().fmtstr + ", bonds={bonds!r}"

    def __eq__(self, other):
        return np.allclose(self.bonds, other.bonds) and super().__eq__(other)

    def __lt__(self, other):
        return (self.bonds < other.bonds and super().__le__(other)) or \
            (self.bonds <= other.bonds and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('bonds')
        return attrs

    @property
    def bonds(self):
        """Return `BondAtom` coordination number."""
        return self._bonds

    @bonds.setter
    def bonds(self, value):
        """Set `BondAtom` bonds."""
        self._bonds = Bonds(value)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(bonds=self.bonds))
        return super_dict


class BondAtoms(Atoms):
    """An `Atoms` sub-class for `BondAtom`\ s.

    A container class for :class:`~sknano.core.atoms.BondAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `BondAtoms`}, optional
        if not `None`, then a list of `BondAtom` instance objects or an
        existing `BondAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return BondAtom

    def sort(self, key=attrgetter('bonds'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def bonds(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`BondAtom.bonds`\ s."""
        return np.asarray([atom.bonds for atom in self])
