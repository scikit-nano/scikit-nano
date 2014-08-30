# -*- coding: utf-8 -*-
"""
===============================================================================
Extended Atoms class feature set (:mod:`sknano.core.atoms._extended_atoms`)
===============================================================================

An "eXtended" `Atoms` class for structure analysis.

.. currentmodule:: sknano.core.atoms._extended_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from functools import total_ordering
from operator import attrgetter

import numpy as np
from sknano.core import xyz
from ._atoms import Atoms

__all__ = ['XAtoms']


@total_ordering
class XAtoms(Atoms):
    """An eXtended `Atoms` class with more atom attributes.

    Parameters
    ----------
    atoms : {None, sequence, `XAtoms`}, optional
        if not `None`, then a list of `XAtom` instance objects or an
        existing `XAtoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False):

        super(XAtoms, self).__init__(atoms=atoms,
                                     copylist=copylist,
                                     deepcopy=deepcopy)
        self._atomtypes = {}

    def __eq__(self, other):
        return self[:] == other[:]

    def __lt__(self, other):
        return self[:] < other[:]

    def sort(self, key=None, reverse=False):
        if key is None:
            self._data.sort(key=attrgetter('element', 'Z', 'atomtype',
                                           'moleculeID', 'atomID'),
                            reverse=reverse)
        else:
            self._data.sort(key=key, reverse=reverse)

    @property
    def atomtypes(self):
        """Return the atom type dictionary."""
        self._update_atomtypes()
        return self._atomtypes

    def _update_atomtypes(self):
        for atom in self:
            if atom.atomtype not in self._atomtypes:
                self._atomtypes[atom.atomtype] = {}
                self._atomtypes[atom.atomtype]['mass'] = atom.m
                self._atomtypes[atom.atomtype]['q'] = atom.q

    @property
    def atom_ids(self):
        """Return array of `XAtom` IDs."""
        return np.asarray([atom.atomID for atom in self])

    @property
    def charges(self):
        """Return array of `XAtom` charges."""
        return np.asarray([atom.q for atom in self])

    @property
    def Ntypes(self):
        """Number of atom types in `XAtoms`."""
        return len(self.atomtypes.keys())

    @property
    def q(self):
        """Return the total net charge of `XAtoms`."""
        return self.charges.sum()

    @property
    def velocities(self):
        """Return array of `XAtom` velocities."""
        return np.asarray([atom.v for atom in self])

    def add_atomtype(self, atom):
        """Add atom type to :attr:`~XAtoms.atomtypes`.

        Parameters
        ----------
        atom : :class:`~sknano.core.atoms.XAtom`
            A :class:`~sknano.core.atoms.XAtom` instance.

        """
        if atom.atomtype not in self._atomtypes:
            self._atomtypes[atom.atomtype] = {}
            self._atomtypes[atom.atomtype]['mass'] = atom.m
            self._atomtypes[atom.atomtype]['q'] = atom.q

    def add_atomtypes(self, atoms=None):
        """Add atomtype for each atom in atoms to atomtypes dictionary.
        Parameters
        ----------
        atoms : sequence
            a list of `XAtom` object instances

        """
        try:
            [self.add_atomtype(atom) for atom in atoms]
        except TypeError:
            print('Expected an iterable sequence of `XAtom` objects.')

    def assign_unique_ids(self, starting_id=1):
        """Assign unique ID to each `XAtom` in `XAtoms`."""
        for i, atom in enumerate(self, start=starting_id):
            atom.atomID = i

    def filter_atoms(self, filter_atom_ids, invert=False):
        """Filter `XAtoms`.

        Parameters
        ----------
        filter_atom_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `XAtoms`

        """
        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        return XAtoms(np.asarray(self)[filter_indices].tolist())

    def get_atom(self, atomID=None, index=None):
        try:
            return self[atomID - 1]
        except (TypeError, IndexError):
            try:
                return self[index]
            except (TypeError, IndexError):
                return None

    def get_filtered_atom_ids(self, filter_atom_ids, invert=False):
        """Return atom ids filtered by list of `filter_atom_ids`.

        Parameters
        ----------
        invert : bool, optional

        Returns
        -------
        ndarray

        """
        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        return self.atom_ids[filter_indices]

    def get_filtered_coords(self, filter_atom_ids, components=None,
                            asdict=False, invert=False):
        """Return filtered coordinates filtered by filter_atom_ids.

        Parameters
        ----------
        filter_atom_ids : array_like
        components : {None, sequence}, optional
        asdict : bool, optional
        invert : bool, optional

        Returns
        -------
        filtered_coords : :py:class:`python:~collections.OrderedDict` or
        ndarray

        """
        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        filtered_coords = self.coords[filter_indices]

        if components is None or components == 'r':
            components = ('x', 'y', 'z')
        elif isinstance(components, (str, unicode)):
            components = (components,)

        if asdict:
            return OrderedDict(zip(components,
                                   [filtered_coords[:, xyz.index(component)]
                                    for component in components]))
        else:
            filtered_coords

    def select(self, **kwargs):
        pass

    def select_within(self, volume):
        pass
