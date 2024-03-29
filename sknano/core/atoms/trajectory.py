# -*- coding: utf-8 -*-
"""
===============================================================================
Trajectory class for MD simulations (:mod:`sknano.core.atoms.trajectory`)
===============================================================================

Classes for analyzing the atom trajectories of molecular dynamics simulations.

.. currentmodule:: sknano.core.atoms.trajectory

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import BaseClass, UserList, TabulateMixin
from .md_atoms import MDAtom as Atom, MDAtoms as Atoms

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
                # snapshot.nselected = snapshot.Natoms
        else:
            snapshot = self.traj.get_snapshot(ts)
            for i in range(snapshot.Natoms):
                snapshot.atom_selection[i] = True
            # snapshot.nselected = snapshot.Natoms

    def update(self, atoms, ts=None):
        atom_ids = atoms.ids
        if ts is None:
            for ss in self.traj:
                if not ss.selected:
                    continue
                aselection = ss.atom_selection
                id_idx = ss.atomattrs.index('id')
                for i, atom in enumerate(ss.get_atoms(asarray=True)):
                    if int(atom[id_idx]) in atom_ids:
                        aselection[i] = True
                    else:
                        aselection[i] = False
        else:
            ss = self.traj.get_snapshot(ts)
            aselection = ss.atom_selection
            id_idx = ss.atomattrs.index('id')
            for i, atom in enumerate(ss.get_atoms(asarray=True)):
                if int(atom[id_idx]) in atom_ids:
                    aselection[i] = True
                else:
                    aselection[i] = False


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
        self.traj.atom_selection.all()
        self.print_fraction_selected()

    def one(self, ts):
        """Select only timestep `ts`."""
        [setattr(snapshot, 'selected', False) for snapshot in self.traj]
        try:
            self.traj.get_snapshot(ts).selected = True
        except AttributeError:
            pass
        self.traj.atom_selection.all()
        self.print_fraction_selected()

    def none(self):
        """Deselect all timesteps."""
        [setattr(snapshot, 'selected', False) for snapshot in self.traj]
        self.print_fraction_selected()

    def skip(self, n):
        """Select every `n`\ th timestep from currently selected timesteps."""
        count = n - 1
        for ss in self.traj:
            if not ss.selected:
                continue
            count += 1
            if count == n:
                count = 0
                continue
            ss.selected = False
        self.traj.atom_selection.all()
        self.print_fraction_selected()

    def print_fraction_selected(self):
        print('{}/{} snapshots selected'.format(
            self.traj.nselected, self.traj.Nsnaps))


class Snapshot(TabulateMixin, BaseClass):
    """Container class for :class:`Trajectory` data at single timestep

    Parameters
    ----------
    trajectory : :class:`Trajectory`, optional

    """
    def __init__(self, trajectory=None):
        super().__init__()

        self.trajectory = trajectory
        self.timestep = None
        self.domain = None
        self._atoms = None
        self._atoms_array = None
        self._formatter = None

        self.fmtstr = "trajectory={trajectory!r}"

    def __getattr__(self, name):
        try:
            return getattr(self.trajectory, name)
        except AttributeError:
            try:
                return getattr(self.formatter, name)
            except AttributeError:
                return super().__getattr__(name)

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        if self.trajectory is not None:
            items = ['timestep', 'Natoms']
            values = [getattr(self, item) for item in items]
            table = self._tabulate(list(zip(items, values)))
            strrep = '\n'.join((strrep, objstr, table))
        return strrep

    @property
    def formatter(self):
        """An instance of :class:`~sknano.io.DUMPFormatter`."""
        return self._formatter

    @formatter.setter
    def formatter(self, value):
        self._formatter = value
        self.atomattrs = value.atomattrs
        self.atomattrmap = value.atomattrmap
        self.attr_dtypes = value.attr_dtypes

    @property
    def atoms(self):
        """Snapshot atoms."""
        if self._atoms is None:
            self._update_atoms()
            return self.atoms
        else:
            return self._atoms.filtered(self.atom_selection)

    @atoms.setter
    def atoms(self, value):
        if isinstance(value, np.ndarray):
            self._atoms_array = value

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
        aselection = self.atom_selection
        return len(aselection[aselection])

    @property
    def nselect(self):
        """Alias for :attr:`Snapshot.nselected`."""
        return self.nselected

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
            return self._atoms_array
        return self.atoms

    def _update_atoms(self):
        atoms = Atoms()
        traj = self.trajectory
        atomattrs = self.atomattrs
        atomattrmap = self.atomattrmap
        attr_dtypes = self.attr_dtypes
        id_idx = atomattrs.index('id')
        for atom in self._atoms_array:
            try:
                reference_atom = \
                    traj.reference_atoms.get_atom(int(atom[id_idx]))
            except AttributeError:
                reference_atom = None

            attrs = {attr: attr_dtypes[idx](atom[idx]) for idx, attr in
                     enumerate(atomattrs)}

            atoms.append(Atom(reference_atom=reference_atom, **attrs))

        if atomattrmap is not None:
            for (from_attr, to_attr), attrmap in atomattrmap.items():
                atoms.mapatomattr(from_attr=from_attr, to_attr=to_attr,
                                  attrmap=attrmap)
        self._atoms = atoms

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(trajectory=self.trajectory)


class Trajectory(TabulateMixin, UserList):
    """Base class for trajectory analysis.

    Parameters
    ----------
    snapshots : :class:`~python:list`, optional

    """
    def __init__(self, snapshots=None):
        super().__init__(initlist=snapshots)
        self.fmtstr = "snapshots={snapshots!r}"
        self.time_selection = TimeSelection(self)
        self.atom_selection = AtomSelection(self)
        self.reference_atoms = None
        self._reference_snapshot = None

    @property
    def __item_class__(self):
        return Snapshot

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        if self.data:
            items = ['Nsnaps', 'nselected', 'timesteps']
            values = [self.Nsnaps, self.nselected, self.timesteps]
            table = self._tabulate(list(zip(items, values)))
            strrep = '\n'.join((strrep, objstr, table))
            if self._reference_snapshot is not None:
                strrep = '\n'.join((strrep, str(self._reference_snapshot)))
        return strrep

    @property
    def Nsnaps(self):
        """Number of :class:`Snapshot`\ s in `Trajectory`."""
        return len(self.data)

    @property
    def atom_selection(self):
        """:class:`AtomSelection` class."""
        return self._atom_selection

    @atom_selection.setter
    def atom_selection(self, value):
        if not isinstance(value, AtomSelection):
            raise ValueError('Expected an `AtomSelection` instance.')
        self._atom_selection = value

    @property
    def time_selection(self):
        """:class:`TimeSelection` class."""
        return self._time_selection

    @time_selection.setter
    def time_selection(self, value):
        if not isinstance(value, TimeSelection):
            raise ValueError('Expected a `TimeSelection instance.')
        self._time_selection = value

    @property
    def aselect(self):
        """Alias for :attr:`~Trajectory.atom_selection`."""
        return self.atom_selection

    @aselect.setter
    def aselect(self, value):
        self.atom_selection = value

    @property
    def tselect(self):
        """Alias for :attr:`~Trajectory.time_selection`."""
        return self.time_selection

    @tselect.setter
    def tselect(self, value):
        self.time_selection = value

    @property
    def nselected(self):
        """Number of selected snapshots."""
        n = 0
        for ss in self:
            if ss.selected:
                n += 1
        return n

    @property
    def nselect(self):
        """Alias for :attr:`~Trajectory.nselected`."""
        return self.nselected

    @property
    def snapshots(self):
        """Returns the list of :class:`Trajectory` :class:`Snapshot`\ s."""
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
        return None

    def timestep_index(self, ts):
        """Return index of :class:`Snapshot` with timestep `ts`."""
        try:
            return self.timesteps.index(ts)
        except ValueError:
            print("No timestep {:d} exists".format(ts))
            return None

    @property
    def reference_snapshot(self):
        """Reference snapshot for computing changes in atom trajectories."""
        return self._reference_snapshot

    @reference_snapshot.setter
    def reference_snapshot(self, value):
        if not isinstance(value, Snapshot):
            raise TypeError('Expected a `Snapshot` instance.')
        self._reference_snapshot = value
        self.reference_atoms = self.reference_snapshot.atoms
        self.reference_atoms.update_attrs()

    @property
    def timesteps(self):
        """List of selected :class:`Trajectory` :class:`Snapshot` \
            :attr:`~Snapshot.timestep`\ s."""
        v = np.zeros(self.nselected, dtype=int)
        for i, snapshot in enumerate(self.data):
            if snapshot.selected:
                v[i] = snapshot.timestep
        return v.tolist()

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(snapshots=self.data)
