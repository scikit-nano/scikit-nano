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

from operator import attrgetter

import numpy as np
from ._atoms import Atoms

__all__ = ['XAtoms']


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
    atomattrs = Atoms.atomattrs + \
        ['atomID', 'moleculeID', 'atomtype', 'q', 'v', 'vx', 'vy', 'vz',
         'n', 'nx', 'ny', 'nz']

    def __init__(self, atoms=None, copylist=True, deepcopy=False):

        super(XAtoms, self).__init__(atoms=atoms,
                                     copylist=copylist,
                                     deepcopy=deepcopy)
        self._atomtypes = {}

    def sort(self, key=None, reverse=False):
        if key is None:
            self.data.sort(key=attrgetter('element', 'Z', 'atomtype',
                                          'moleculeID', 'atomID'),
                           reverse=reverse)
        else:
            self.data.sort(key=key, reverse=reverse)

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
        if len(set([atom.atomID for atom in self])) != self.Natoms:
            self.assign_unique_ids()
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

    def filter_ids(self, atom_ids, invert=False):
        """Return `Atoms` by :attr:`XAtoms.atom_ids` in `atom_ids`.

        Parameters
        ----------
        filter_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `XAtoms`

        """
        filtered_atoms = \
            np.asarray(self)[np.in1d(
                self.atom_ids, atom_ids, invert=invert).nonzero()].tolist()
        return self.__class__(atoms=filtered_atoms)

    def get_atom(self, atomID):
        try:
            return self[np.where(self.atom_ids == atomID)[0]]
        except TypeError:
            print('No atom with atomID = {}'.format(atomID))
            return None

    def getatomattr(self, attr):
        if attr not in self.atomattrs:
            errmsg = '{} not in list of allowed attributes:\n{}'
            raise ValueError(errmsg.format(attr, self.atomattrs))
        return np.asarray([getattr(atom, attr) for atom in self])

    #def select(self, cmd):
    #    cmdsplit =
    #    if cmd.lower().startswith('with attr'):
    #        if attr not in self.attributes:
    #            errmsg = '{} not in list of allowed attributes:\n{}'
    #            raise ValueError(errmsg.format(attr, self.attributes))
    #        return np.asarray([getattr(atom, attr) for atom in self])

    def select_within(self, volume):
        pass
