# -*- coding: utf-8 -*-
"""
===============================================================================
Trajectory class for MD simulations (:mod:`sknano.core.atoms._trajectory`)
===============================================================================

Classes for analyzing the atom trajectories of molecular dynamics simulations.

.. currentmodule:: sknano.core.atoms._trajectory

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import BaseClass, UserList
from ._md_atoms import MDAtom as Atom, MDAtoms as Atoms

__all__ = ['Snapshot', 'Trajectory']


class AtomSelection:
    """:class:`Trajectory` atom selection class.

    Parameters
    ----------
    traj : :class:`Trajectory`

    """
    def __init__(self, traj):
        self.traj = traj

    def all(self, ts=None):
        """Select all atoms for all snapshots or snapshot at given timestep.

        Parameters
        ----------
        ts : {None, int}, optional

        """
        if ts is None:
            for snapshot in self.traj:
                if not snapshot.selected:
                    continue
                for i in range(snapshot.Natoms):
                    snapshot.atom_selection[i] = True
                snapshot.nselected = snapshot.Natoms
        else:
            snapshot = self.traj.get_snapshot(ts)
            for i in range(snapshot.Natoms):
                snapshot.atom_selection[i] = True
            snapshot.nselected = snapshot.Natoms


class TimeSelection:
    """:class:`Trajectory` time selection class.

    Parameters
    ----------
    traj : :class:`Trajectory`

    """
    def __init__(self, traj):
        self.traj = traj

    def all(self, ts=None):
        """Select all trajectory snapshots/timesteps."""
        [setattr(snapshot, 'selected', True) for snapshot in self.traj]
        self.traj.nselected = self.traj.Nsnaps
        self.traj.atom_selection.all()
        self.print_fraction_selected()

    def one(self, ts):
        """Select only timestep `ts`."""
        [setattr(snapshot, 'selected', False) for snapshot in self.traj]
        try:
            self.traj.get_snapshot(ts).selected = True
            self.traj.nselected = 1
        except AttributeError:
            pass
        self.traj.atom_selection.all()
        self.print_fraction_selected()

    def none(self):
        """Deselect all timesteps."""
        [setattr(snapshot, 'selected', False) for snapshot in self.traj]
        self.traj.nselected = 0
        self.print_fraction_selected()

    def skip(self, n):
        """Select every `n`\ th timestep from currently selected timesteps."""
        count = n - 1
        for snapshot in self.traj:
            if not snapshot.selected:
                continue
            count += 1
            if count == n:
                count = 0
                continue
            snapshot.selected = False
            self.traj.nselected -= 1
        self.traj.atom_selection.all()
        self.print_fraction_selected()

    def print_fraction_selected(self):
        print('{}/{} snapshots selected'.format(
            self.traj.nselected, self.traj.Nsnaps))


class Snapshot(BaseClass):
    """Container class for :class:`Trajectory` data at single timestep"""
    def __init__(self, trajectory=None):

        super().__init__()

        self.trajectory = trajectory

        self.atomattrs = None
        self.attr_dtypes = None
        self.timestep = None

        self._atoms = None

        self.fmtstr = "trajectory={trajectory!r}"

    @property
    def atoms(self):
        """Snapshot atoms."""
        atoms = Atoms()
        for atom in self._atoms:
            try:
                reference_atom = \
                    self.trajectory.reference_atoms.get_atom(
                        int(atom[self.atomattrs.index('id')]))
            except AttributeError:
                reference_atom = None

            try:
                t0_atom = self.trajectory.t0_atoms.get_atom(
                    int(atom[self.atomattrs.index('id')]))
            except AttributeError:
                t0_atom = None

            attrs = [dtype(value) for dtype, value in
                     zip(self.attr_dtypes, atom)]
            atoms.append(Atom(reference_atom=reference_atom,
                              t0_atom=t0_atom,
                              **dict(list(zip(self.atomattrs, attrs)))))
        return atoms

    @atoms.setter
    def atoms(self, value):
        self._atoms = value

    @property
    def atom_selection(self):
        """:class:`~numpy:numpy.ndarray` boolean array."""
        return self._atom_selection

    @atom_selection.setter
    def atom_selection(self, value):
        if not isinstance(value, (list, np.ndarray)):
            raise ValueError('Expected an array_like object.')
        self._atom_selection = np.asarray(value, dtype=bool)

    @property
    def aselect(self):
        """Alias for :attr:`Snapshot.atom_selection`."""
        return self.atom_selection

    @aselect.setter
    def aselect(self, value):
        self.atom_selection = value

    @property
    def selected(self):
        """True/False if this snapshot is selected."""
        return self._selected

    @selected.setter
    def selected(self, value):
        self._selected = bool(value)

    @property
    def tselect(self):
        """Alias for :attr:`Snapshot.selected`."""
        return self.selected

    @tselect.setter
    def tselect(self, value):
        self.selected = value

    @property
    def nselected(self):
        """Number of selected atoms in this snapshot."""
        return self._nselected

    @nselected.setter
    def nselected(self, value):
        self._nselected = int(value)

    @property
    def nselect(self):
        """Alias for :attr:`Snapshot.nselected`."""
        return self.nselected

    @nselect.setter
    def nselect(self, value):
        self.nselected = value

    def get_atoms(self, asarray=False):
        """Get atoms.

        Parameters
        ----------
        asarray : :class:`~python:bool`

        Returns
        -------
        :class:`~numpy:numpy.ndarray` or :class:`MDAtoms`
            if `asarray` is `True`, the atoms are returned as an
            :class:`~numpy:numpy.ndarray`, otherwise an :class:`MDAtoms`
            instance is returned.

        """
        if asarray:
            return self._atoms
        return self.atoms

    def todict(self):
        return dict(trajectory=self.trajectory)


class Trajectory(BaseClass, UserList):
    """Base class for trajectory analysis."""

    def __init__(self, snapshots=None):
        super().__init__(initlist=snapshots)
        self.fmtstr = "snapshots={snapshots!r}"
        self.time_selection = TimeSelection(self)
        self.atom_selection = AtomSelection(self)
        self.nselected = 0
        self.reference_atoms = None
        self._reference_snapshot = None

        self.t0_atoms = None
        self._t0_snapshot = None

    @property
    def Nsnaps(self):
        """Number of :class:`Snapshot`\ s in `Trajectory`."""
        return len(self.data)

    @property
    def atom_selection(self):
        """`AtomSelection` class."""
        return self._atom_selection

    @atom_selection.setter
    def atom_selection(self, value):
        if not isinstance(value, AtomSelection):
            raise ValueError('Expected an `AtomSelection` instance.')
        self._atom_selection = value

    @property
    def time_selection(self):
        return self._time_selection

    @time_selection.setter
    def time_selection(self, value):
        if not isinstance(value, TimeSelection):
            raise ValueError('Expected a `TimeSelection instance.')
        self._time_selection = value

    @property
    def aselect(self):
        """Alias for :attr:`Trajectory.atom_selection`."""
        return self.atom_selection

    @aselect.setter
    def aselect(self, value):
        self.atom_selection = value

    @property
    def tselect(self):
        """Alias for :attr:`Trajectory.time_selection`."""
        return self.time_selection

    @tselect.setter
    def tselect(self, value):
        self.time_selection = value

    @property
    def nselected(self):
        """Number of selected snapshots."""
        return self._nselected

    @nselected.setter
    def nselected(self, value):
        self._nselected = int(value)

    @property
    def nselect(self):
        """Alias for :attr:`Trajectory.nselected`."""
        return self.nselected

    @nselect.setter
    def nselect(self, value):
        self.nselected = value

    @property
    def snapshots(self):
        """Returns the list of :class:`Snapshot`\ s."""
        return self.data

    def sort(self, key=attrgetter('timestep'), reverse=False):
        """Sort the trajectory :class:`Snapshot`\ s."""
        super().sort(key=key, reverse=reverse)

    def cull(self):
        """Remove duplicate timesteps from `Trajectory`."""
        i = 1
        while i < len(self.data):
            if self.data[i].timestep == self.data[i-1].timestep:
                del self.data[i]
            else:
                i += 1

    def get_snapshot(self, ts):
        """Return :class:`Snapshot` with timestep `ts`."""
        for snapshot in self:
            if snapshot.timestep == ts:
                return snapshot
        print("No snapshot at ts={:d} exists".format(ts))

    def timestep_index(self, ts):
        """Return index of :class:`Snapshot` with timestep `ts`."""
        for i, snapshot in enumerate(self):
            if snapshot.timestep == ts:
                return i
        print("No timestep {:d} exists".format(ts))

    @property
    def reference_snapshot(self):
        return self._reference_snapshot

    @reference_snapshot.setter
    def reference_snapshot(self, value):
        if not isinstance(value, Snapshot):
            raise TypeError('Expected a `Snapshot` instance.')
        self._reference_snapshot = value
        self.reference_atoms = self.reference_snapshot.atoms
        self.reference_atoms.update_attrs()

    @property
    def t0_snapshot(self):
        return self._t0_snapshot

    @t0_snapshot.setter
    def t0_snapshot(self, value):
        if not isinstance(value, Snapshot):
            raise TypeError('Expected a `Snapshot` instance.')
        self._t0_snapshot = value
        self.t0_atoms = self.t0_snapshot.atoms
        self.t0_atoms.update_attrs()

    @property
    def timesteps(self):
        v = np.zeros(self.nselected, dtype=int)
        for i, snapshot in enumerate(self.data):
            if snapshot.selected:
                v[i] = snapshot.timestep
        return v

    def todict(self):
        return dict(snapshots=self.data)
