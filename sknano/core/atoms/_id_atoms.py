# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with id attributes (:mod:`sknano.core.atoms._id_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._id_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter
import numbers

import numpy as np

from ._atoms import Atom, Atoms

__all__ = ['IDAtom', 'IDAtoms']


@total_ordering
class IDAtom(Atom):
    """An `Atom` class with id attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    id : int, optional
        atom ID
    mol : int, optional
        molecule ID
    """
    def __init__(self, *args, id=0, mol=0, **kwargs):

        if 'atomID' in kwargs:
            id = kwargs['atomID']
            del kwargs['atomID']

        if 'moleculeID' in kwargs:
            mol = kwargs['moleculeID']
            del kwargs['moleculeID']

        super().__init__(*args, **kwargs)

        self.id = id
        self.mol = mol
        self.fmtstr = super().fmtstr + ", id={id!r}, mol={mol!r}"

    def __eq__(self, other):
        return self.id == other.id and self.mol == other.mol and \
            super().__eq__(other)

    def __lt__(self, other):
        return ((self.id < other.id and self.mol <= other.mol and
                 super().__le__(other))
                or (self.id <= other.id and self.mol < other.mol and
                    super().__le__(other))
                or (self.id <= other.id and self.mol <= other.mol and
                    super().__lt__(other)))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['id', 'mol'])
        return attrs

    @property
    def id(self):
        """:attr:`~IDAtom.id`."""
        return self._id

    @id.setter
    def id(self, value):
        """Set atom id.

        Parameters
        ----------
        value : int
            atom ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._id = int(value)

    @property
    def mol(self):
        """:attr:`~IDAtom.mol`."""
        return self._mol

    @mol.setter
    def mol(self, value):
        """Set :attr:`~IDAtom.mol`.

        Parameters
        ----------
        value : int
            molecule ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._mol = int(value)

    @property
    def atomID(self):
        return self.id

    @atomID.setter
    def atomID(self, value):
        self.id = value

    @property
    def moleculeID(self):
        return self.mol

    @moleculeID.setter
    def moleculeID(self, value):
        self.mol = value

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(id=self.id, mol=self.mol))
        return super_dict


class IDAtoms(Atoms):
    """An `Atoms` sub-class for `IDAtom`\ s.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.IDAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `IDAtoms`}, optional
        if not `None`, then a list of `IDAtom` instance objects or an
        existing `IDAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return IDAtom

    def sort(self, key=attrgetter('mol', 'id'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def ids(self):
        """Return array of `IDAtom.id`\ s."""
        if len(set([atom.id for atom in self])) != self.Natoms:
            self.assign_unique_ids()
        return np.asarray([atom.id for atom in self])

    @property
    def mols(self):
        """Return array of `IDAtom.mol`\ s."""
        return np.asarray([atom.mol for atom in self])

    @property
    def atom_ids(self):
        """Alias for :attr:`IDAtoms.ids`."""
        return self.ids

    @property
    def mol_ids(self):
        """Alias for :attr:`IDAtoms.mols`."""
        return self.mols

    @property
    def molecule_ids(self):
        """Alias for :attr:`IDAtoms.mols`."""
        return self.mols

    def assign_unique_ids(self, starting_id=1):
        """Assign unique :attr:`IDAtom.id` to each `IDAtom` in `IDAtoms`."""
        [setattr(atom, 'id', i) for i, atom in
         enumerate(self, start=starting_id)]

    def filter_ids(self, atom_ids, invert=False):
        """Filter `Atoms` by :attr:`IDAtoms.ids` in `atom_ids`.

        .. versionchanged:: 0.3.11

           Filters `Atoms` **in-place**. Use :meth:`~IDAtoms.filtered_ids`
           to get a new list of `Atoms`.

        Parameters
        ----------
        atom_ids : array_like
        invert : bool, optional

        """
        mask = np.in1d(self.ids, atom_ids, invert=invert).nonzero()
        self.data = np.asarray(self)[mask].tolist()

    def filtered_ids(self, atom_ids, invert=False):
        """Return new `Atoms` filtered by :attr:`IDAtoms.ids` in `atom_ids`.

        .. versionadded:: 0.3.11

        Parameters
        ----------
        atom_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `Atoms`
            An instance of `Atoms` (sub)class.

        """
        mask = np.in1d(self.ids, atom_ids, invert=invert).nonzero()
        return self.__class__(atoms=np.asarray(self)[mask].tolist(),
                              **self.kwargs)

    def get_atom(self, id):
        """Get `IDAtom` with :attr:`Xatom.id` == `id`.

        Parameters
        ----------
        id : int

        Returns
        -------
        atom : `IDAtom` or `None`
            `IDAtom` instance if `IDAtoms` contains `IDAtom` with
            :attr:`IDAtom.id` == `id`, otherwise `None`

        """
        try:
            return self[np.where(self.ids == id)[0]]
        except TypeError:
            print('No atom with id = {}'.format(id))
            return None
