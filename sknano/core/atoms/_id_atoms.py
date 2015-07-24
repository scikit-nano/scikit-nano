# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for `IDAtom`\ s (:mod:`sknano.core.atoms._id_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._id_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from ._atoms import Atoms
from ._id_atom import IDAtom

__all__ = ['IDAtoms']


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
