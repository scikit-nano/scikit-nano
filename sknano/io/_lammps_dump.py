# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS dump format (:mod:`sknano.io._lammps_dump`)
====================================================================

.. currentmodule:: sknano.io._lammps_dump

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from glob import glob
from operator import attrgetter
# import re
import sys

import numpy as np

from monty.io import zopen
from sknano.core import deprecate_kwarg, get_fpath, flatten, grouper
from sknano.core.atoms import Trajectory, Snapshot, Atoms, MDAtoms, \
    MDAtom as Atom
# from sknano.core.crystallography import Crystal3DLattice
from ._base import StructureData, StructureDataError, StructureDataFormatter

__all__ = ['DUMP', 'DUMPData', 'DUMPReader', 'DUMPWriter', 'DUMPError',
           'DUMPFormatter', 'DUMPIO', 'DUMPIOReader', 'DUMPIOWriter',
           'DUMPIOError', 'DUMPIOFormatter']


class DUMPReader(StructureData):
    """Class for reading `LAMMPS dump` file format.

    Parameters
    ----------
    *args : :class:`~python:list`
        :class:`~python:list` of one or more LAMMPS dump files.
    autoread : :class:`~python:bool`, optional
        Automatically read dump files. Default is `True`

        .. versionadded:: 0.4.0

    reference_timestep : :class:`~python:int`, optional
        The `timestep` corresponding to the
        :attr:`~sknano.core.atoms.Snapshot.timestep` of the
        :class:`~sknano.core.atoms.Trajectory`\ s
        :attr:`~sknano.core.atoms.Trajectory.reference_snapshot`.
        Default is `None`.

        .. note::
           **Overrides** `reference_index` if not `None`

        .. versionadded:: 0.4.0

    reference_index : :class:`~python:int`, optional
        The :class:`~python:list` index corresponding to the
        :attr:`~sknano.core.atoms.Snapshot.timestep` of the
        :class:`~sknano.core.atoms.Trajectory`\ s
        :attr:`~sknano.core.atoms.Trajectory.reference_snapshot`.
        Default is `None`.

        .. versionadded:: 0.4.0

    dumpattrmap : class:`~python:dict`
        Python :class:`~python:dict` mapping custom dump attributes
        to :class:`~sknano.core.atoms.Atom` attributes.
    atomattrmap : :class:`~python:dict`
        :class:`~python:dict` mapping atom type to atom element

    Attributes
    ----------
    trajectory

    Examples
    --------
    Using LAMMPS, one often saves molecular dynamics atom trajectories
    to one or more custom dump files containing per-atom attributes of
    of each atom at a given timestep. Often the per-atom attributes
    are calculated from custom LAMMPS compute commands that one may
    want to map to a specific :class:`~sknano.core.atoms.Atom`
    attribute.

    For example, consider a dump file containing per-atom dump attributes
    *id*, *mol*, *type*, *x*, *y*, *z*, *vx*, *vy*, *vz*, *c_atom_ke*, and
    *c_atom_pe*, where *c_atom_ke* and *c_atom_pe* are the LAMMPS
    *compute-IDs* from custom compute commands defined in a LAMMPS input
    script, which calculate per-atom values for kinetic energy and
    potential energy, respectively. To map these compute-IDs and their
    per-atom values to their respective :class:`~sknano.core.atoms.Atom`
    attributes :attr:`~sknano.core.atom.EnergyAtom.ke` and
    :attr:`~sknano.core.atom.EnergyAtom.pe`, one
    would use the following :class:`~python:dict` mapping::

    >>> dumpattrmap = {'c_atom_pe': 'pe', 'c_atom_ke': 'ke'}

    >>> from sknano.io import DUMPReader
    >>> dumps = DUMPReader('dump.*', dumpattrmap=dumpattrmap)
    >>> print(repr(dumps.dumpattrs2str()))
    'id mol type x y z vx vy vz c_atom_ke c_atom_pe c_atom_CN'
    >>> print(repr(dumps.atomattrs2str()))
    'id mol type x y z vx vy vz ke pe'

    Futhermore, atom attributes such as mass or element symbol are not
    typically included in a LAMMPS dump. The LAMMPS attribute usually
    associated with a specic atom element/mass is the `type` attribute.
    The `atomattrmap` parameter takes a :class:`~python:dict`, mapping
    keys of 2-tuples of strings of the form (`from_attr`, `to_attr`)
    to :class:`~python:dict` values of the form
    {`from_attr_value_1`: `to_attr_value_1`, ...,
    `from_attr_value_N`: `to_attr_value_N`}, where
    `(from_attr|to_attr)_value_(i...N)` are the explicit key, value
    pairs for the (`from_attr`, `to_attr`) attribute names.

    For example, suppose that we want to use the dump attribute `type` values
    to set the :attr:`~sknano.core.atoms.Atom.element` values, where we
    have knowledge of the fact that atoms with `type=1` correspond to
    Carbon atoms (i.e. `element='C'`), while atoms with `type=2` correspond to
    Nitrogen atoms (i.e., `element='N'). We would like to pass this metadata
    to the `DUMPReader` constructor so that the
    :class:`~sknano.core.atoms.Atoms` object returned by the
    :class:`~sknano.core.atoms.Trajectory`
    :class:`~sknano.core.atoms.Snapshot`\ s
    :attr:`~sknano.core.atoms.Snapshot.atoms` attribute will have
    the correct :attr:`~sknano.core.atoms.Atoms.elements` corresponding to the
    atom :attr:`~sknano.core.atoms.TypeAtoms.types`. To achieve this,
    we would call :class:`DUMPReader` with the `atomattrmap` keyword
    argument set as follow::

    >>> dumps = DUMPReader('dump.*', dumpattrmap=dumpattrmap,
    ...                    atomattrmap={('type', 'element'): {1:'C', 2:'Ar'}})
    >>> atoms = dumps[0].atoms
    >>> print(atoms[0])

    """
    @deprecate_kwarg(kwarg='attrmap', since='0.4.0', alternative='dumpattrmap')
    def __init__(self, *args, autoread=True, reference_timestep=None,
                 reference_index=None, formatter=None, style=None,
                 dumpattrs=None, dumpattrmap=None, atomattrmap=None,
                 **kwargs):

        # if 'attrmap' in kwargs:
        #     msg = ("The {!s} keyword argument `attrmap` was deprecated in "
        #            "in version 0.4.0.\nUse `dumpattrmap` instead."
        #            .format(self.__class__.__name__))
        #     warnings.warn(msg, DeprecationWarning, stacklevel=2)
        #     dumpattrmap = kwargs.pop('attrmap')

        super().__init__(**kwargs)
        self.trajectory = Trajectory()
        self.dumpfiles = tuple(flatten([[glob(f) for f in arg.split()]
                                        for arg in args]))
        self._reference_timestep = reference_timestep
        self._reference_index = reference_index

        if formatter is None or not isinstance(formatter, DUMPFormatter):
            formatter = DUMPFormatter(style=style,
                                      dumpattrs=dumpattrs,
                                      dumpattrmap=dumpattrmap,
                                      atomattrmap=atomattrmap)

        self.formatter = formatter

        self.fmtstr = "{dumpfiles!r}, autoread=True, " + \
            "reference_timestep={reference_timestep!r}, " + \
            "reference_index={reference_index!r}, " + formatter.fmtstr

        if autoread and len(self.dumpfiles) > 0:
            self.read()

    def __getattr__(self, name):
        try:
            return getattr(self.trajectory, name)
        except AttributeError:
            try:
                return getattr(self.formatter, name)
            except AttributeError:
                return super().__getattr__(name)

    def __getitem__(self, index):
        return self.trajectory[index]

    def __iter__(self):
        return iter(self.trajectory)

    @property
    def reference_index(self):
        """Reference snapshot index."""
        return self._reference_index

    @reference_index.setter
    def reference_index(self, value):
        self._reference_index = value
        self._update_reference_snapshot()

    def _update_reference_snapshot(self):
        try:
            if len(self.trajectory) > 1:
                if all([value is None for value in
                        (self.reference_timestep, self.reference_index)]):
                    self._reference_index = 0
                    self._update_reference_snapshot()
                elif self.reference_timestep is not None and \
                        self.reference_index is None:
                    index = \
                        self.trajectory.timestep_index(self.reference_timestep)
                    if index is None:
                        index = 0
                        self._reference_timestep = None
                    self._reference_index = index
                    self._update_reference_snapshot()
                elif self.reference_timestep is None and \
                        self.reference_index is not None:
                    try:
                        self._reference_timestep = \
                            self.trajectory.timesteps[self.reference_index]
                    except IndexError:
                        self._reference_index = 0
                    self._update_reference_snapshot()
                else:
                    self.trajectory.reference_snapshot = \
                        self.trajectory[self.reference_index]
        except IndexError:
            pass

    @property
    def reference_timestep(self):
        """Reference snapshot timestep."""
        return self._reference_timestep

    @reference_timestep.setter
    def reference_timestep(self, value):
        self._reference_timestep = value
        self._update_reference_snapshot()

    def read(self):
        """Read all snapshots from each dump file."""
        trajectory = self.trajectory
        for dumpfile in self.dumpfiles:
            with zopen(dumpfile) as f:
                snapshot = self.read_snapshot(f)
                while snapshot is not None:
                    trajectory.append(snapshot)
                    print(snapshot.timestep, end=' ')
                    sys.stdout.flush()
                    snapshot = self.read_snapshot(f)
        print()

        trajectory.sort(key=attrgetter('timestep'))
        trajectory.cull()

        print("read {:d} snapshots".format(self.Nsnaps))

        trajectory.time_selection.all()
        self._update_reference_snapshot()

        fmt = self.formatter
        dumpattrs = fmt.dumpattrs

        if dumpattrs:
            print('Dumped Atom attributes: {}'.format(fmt.dumpattrs2str()))
        else:
            print('No dump column assignments')

        if 'x' not in dumpattrs or 'y' not in dumpattrs or \
                'z' not in dumpattrs:
            print('dump scaling status unknown')
        elif self.Nsnaps > 0:
            if self.scale_original:
                self.unscale()
            elif self.scale_original is None:
                print('dump scaling status unknown')
            else:
                print('dump is already unscaled')

    def read_snapshot(self, f):
        """Read snapshot from file."""
        try:
            snapshot = Snapshot(self.trajectory)
            domain = snapshot.domain

            f.readline()
            snapshot.timestep = int(f.readline().strip().split()[0])
            f.readline()
            Natoms = snapshot.Natoms = int(f.readline().strip())
            snapshot.atom_selection = np.zeros(Natoms, dtype=bool)

            item = f.readline().strip()
            try:
                snapshot.boxstr = item.split('BOUNDS')[1].strip()
            except IndexError:
                snapshot.boxstr = ''

            if 'xy' in snapshot.boxstr:
                domain.triclinic = True

            for dim, tilt_factor in zip(('x', 'y', 'z'), ('xy', 'xz', 'yz')):
                bounds = f.readline().strip().split()

                setattr(domain, dim + 'lo', float(bounds[0]))
                setattr(domain, dim + 'hi', float(bounds[1]))

                if domain.triclinic:
                    setattr(domain, tilt_factor, float(bounds[2]))

            if domain.triclinic:
                xlo_bound = domain.xlo
                xhi_bound = domain.xhi
                ylo_bound = domain.ylo
                yhi_bound = domain.yhi
                xy = domain.xy
                xz = domain.xz
                yz = domain.yz
                domain.xlo = xlo_bound - min((0.0, xy, xz, xy + xz))
                domain.xhi = xhi_bound - max((0.0, xy, xz, xy + xz))
                domain.ylo = ylo_bound - min((0.0, yz))
                domain.yhi = yhi_bound - max((0.0, yz))

            formatter = self.formatter
            if formatter.dumpattrs is None:
                dumpattrs2index = formatter.dumpattrs2index
                xflag = yflag = zflag = None
                attrs = f.readline().strip().split()[2:]
                for i, attr in enumerate(attrs):
                    if attr in ('x', 'xu', 'xs', 'xsu'):
                        dumpattrs2index['x'] = i
                        if attr in ('x', 'xu'):
                            xflag = False
                        else:
                            xflag = True
                    elif attr in ('y', 'yu', 'ys', 'ysu'):
                        dumpattrs2index['y'] = i
                        if attr in ('y', 'yu'):
                            yflag = False
                        else:
                            yflag = True
                    elif attr in ('z', 'zu', 'zs', 'zsu'):
                        dumpattrs2index['z'] = i
                        if attr in ('z', 'zu'):
                            zflag = False
                        else:
                            zflag = True
                    else:
                        dumpattrs2index[attr] = i

                self.scale_original = None
                if all([flag is False for flag in (xflag, yflag, zflag)]):
                    self.scale_original = False
                if all([flag for flag in (xflag, yflag, zflag)]):
                    self.scale_original = True

            else:
                f.readline()

            atoms_array = \
                np.zeros((Natoms, len(formatter.dumpattrs)), dtype=float)
            for n in range(Natoms):
                line = [float(value) for col, value in
                        enumerate(f.readline().strip().split())]
                atoms_array[n] = line

            snapshot._atoms_array = atoms_array
            snapshot.formatter = formatter
            return snapshot
        except IndexError:
            return None

    def scale(self):
        """Scale cartesian coordinates to fractional coordinates."""
        dumpattrs2index = self.formatter.dumpattrs2index
        xi = dumpattrs2index['x']
        yi = dumpattrs2index['y']
        zi = dumpattrs2index['z']
        for snapshot in self.trajectory:
            atoms = snapshot.get_atoms(asarray=True)
            if atoms is not None:
                domain = snapshot.domain
                lx = domain.lx
                ly = domain.ly
                lz = domain.lz
                xy = domain.xy
                xz = domain.xz
                yz = domain.yz

                if np.allclose([xy, xz, yz], np.zeros(3)):
                    atoms[:, xi] = (atoms[:, xi] - domain.xlo) / lx
                    atoms[:, yi] = (atoms[:, yi] - domain.ylo) / ly
                    atoms[:, zi] = (atoms[:, zi] - domain.zlo) / lz
                else:
                    xlo_bound = domain.xlo_bound
                    ylo_bound = domain.ylo_bound
                    zlo_bound = domain.zlo_bound
                    atoms[:, xi] = (atoms[:, xi] - xlo_bound) / lx + \
                        (atoms[:, yi] - ylo_bound) * xy / (lx * ly) + \
                        (atoms[:, zi] - zlo_bound) * (yz * xy - ly * xz) / \
                        (lx * ly * lz)
                    atoms[:, yi] = (atoms[:, yi] - ylo_bound) / ly + \
                        (atoms[:, zi] - zlo_bound) * yz / (ly * lz)
                    atoms[:, zi] = (atoms[:, zi] - zlo_bound) / lz

    def unscale(self):
        """Unscale fractional coordinates to cartesian coordinates."""
        dumpattrs2index = self.formatter.dumpattrs2index
        xi = dumpattrs2index['x']
        yi = dumpattrs2index['y']
        zi = dumpattrs2index['z']
        for snapshot in self.trajectory:
            atoms = snapshot.get_atoms(asarray=True)
            if atoms is not None:
                domain = snapshot.domain
                lx = domain.lx
                ly = domain.ly
                lz = domain.lz
                xy = domain.xy
                xz = domain.xz
                yz = domain.yz
                if np.allclose([xy, xz, yz], np.zeros(3)):
                    atoms[:, xi] = domain.xlo + atoms[:, xi] * lx
                    atoms[:, yi] = domain.ylo + atoms[:, yi] * ly
                    atoms[:, zi] = domain.zlo + atoms[:, zi] * lz
                else:
                    xlo_bound = domain.xlo_bound
                    ylo_bound = domain.ylo_bound
                    zlo_bound = domain.zlo_bound

                    atoms[:, xi] = xlo_bound + atoms[:, xi] * lx + \
                        atoms[:, yi] * xy + atoms[:, zi] * xz
                    atoms[:, yi] = ylo_bound + atoms[:, yi] * ly + \
                        atoms[:, zi] * yz
                    atoms[:, zi] = zlo_bound + atoms[:, zi] * lz

    def wrap(self):
        """Wrap coordinates from outside box to inside."""
        dumpattrs2index = self.formatter.dumpattrs2index
        x = dumpattrs2index['x']
        y = dumpattrs2index['y']
        z = dumpattrs2index['z']

        ix = dumpattrs2index.get('ix', None)
        iy = dumpattrs2index.get('iy', None)
        iz = dumpattrs2index.get('iz', None)

        if all([iflag is not None for iflag in (ix, iy, iz)]):
            for snapshot in self.trajectory:
                atoms = snapshot.get_atoms(asarray=True)
                if atoms is not None:
                    domain = snapshot.domain
                    lx = domain.lx
                    ly = domain.ly
                    lz = domain.lz
                    atoms[:, x] -= atoms[:, ix] * lx
                    atoms[:, y] -= atoms[:, iy] * ly
                    atoms[:, z] -= atoms[:, iz] * lz

    def unwrap(self):
        """Unwrap coordinates from inside box to outside."""
        dumpattrs2index = self.formatter.dumpattrs2index
        x = dumpattrs2index['x']
        y = dumpattrs2index['y']
        z = dumpattrs2index['z']

        ix = dumpattrs2index.get('ix', None)
        iy = dumpattrs2index.get('iy', None)
        iz = dumpattrs2index.get('iz', None)

        if all([iflag is not None for iflag in (ix, iy, iz)]):
            for snapshot in self.trajectory:
                atoms = snapshot.get_atoms(asarray=True)
                if atoms is not None:
                    domain = snapshot.domain
                    lx = domain.lx
                    ly = domain.ly
                    lz = domain.lz
                    atoms[:, x] += atoms[:, ix] * lx
                    atoms[:, y] += atoms[:, iy] * ly
                    atoms[:, z] += atoms[:, iz] * lz

    def new_dumpattr(self, attr, values=None):
        """Add new dump attr to :attr:`~DUMPReader.trajectory.snapshots`.

        This method adds a new dump attribute to the
        :attr:`DUMPFormatter.dumpattrs` :class:`~python:list` and
        adds a new column of data to the
        :attr:`~DUMPReader.trajectory.snapshots._atoms_array` with
        values of 0.0 if `values=None`.

        This method has no effect if `attr` already exists.

        Parameters
        ----------
        attr : :class:`~python:str`
        values : {:class:`~python:float`, array_like}, optional

        """
        dumpattrs2index = self.formatter.dumpattrs2index
        if attr not in dumpattrs2index:
            attridx = len(dumpattrs2index)
            dumpattrs2index[attr] = attridx
            trajectory = self.trajectory
            Nsnaps = trajectory.Nsnaps
            if values is None:
                values = 0.0
            elif isinstance(values, str) and values == 'eval':
                values = {attr: attr}

            if np.isscalar(values):
                values = [ss.Natoms * [values] for ss in trajectory]
            elif isinstance(values, (list, np.ndarray)):
                values = np.asarray(values)
                if values.ndim == 1:
                    values = Nsnaps * [values.tolist()]
            elif isinstance(values, dict):
                dumpattrmap = self.dumpattrmap
                if dumpattrmap is not None:
                    dumpattrmap.update(values)
                values = [[atom.getattr(values[attr], 0.0, recursive=True)
                           for atom in ss.atoms] for ss in trajectory]
            values = np.asarray(values)
            for i, ss in enumerate(trajectory):
                ss._atoms_array = \
                    np.insert(ss._atoms_array, attridx, values[i])

    def map(self, *pairs):
        """Update :attr:`~DUMPFormatter.dumpattrs2index` mapping.

        This method is defined for compatibility with the LAMMPS pizza.py
        dump module.

        Parameters
        ----------
        pairs : :class:`~python:tuple`
            2 :class:`~python:tuple`\s of (column number, column name) pairs.

        Notes
        -----
        This function maps column numbers to dump attribute names,
        **not** the column index.

        """
        pairs = tuple(flatten(pairs))
        if len(pairs) % 2 != 0:
            raise ValueError('Expected (column number, column name) pairs')
        self.dumpattrs2index.update({attr: i - 1 for i, attr in
                                     list(grouper(pairs, 2))})

    def newcolumn(self, name):
        """An alias for :meth:`DUMPReader.new_dumpattr` for compatibility \
            LAMMPS pizza.py dump module."""
        self.new_dumpattr(name)

    def update_dumpattr(self, attr, values=None):
        """Add new dump attr to :attr:`~DUMPReader.trajectory.snapshots`.

        This method adds a new dump attribute to the
        :attr:`DUMPFormatter.dumpattrs` :class:`~python:list` and
        adds a new column of data to the
        :attr:`~DUMPReader.trajectory.snapshots._atoms_array` with
        values of 0.0 if `values=None`.

        This method has no effect if `attr` already exists.

        Parameters
        ----------
        attr : :class:`~python:str`
        values : {:class:`~python:float`, array_like}, optional

        """
        dumpattrs2index = self.formatter.dumpattrs2index
        if attr in dumpattrs2index and values is not None:
            attridx = dumpattrs2index[attr]
            trajectory = self.trajectory
            Nsnaps = trajectory.Nsnaps
            if isinstance(values, str) and values == 'attr':
                values = {attr: attr}
            if np.isscalar(values):
                values = [ss.Natoms * [values] for ss in trajectory]
            elif isinstance(values, (list, np.ndarray)):
                values = np.asarray(values)
                if values.ndim == 1:
                    values = Nsnaps * [values.tolist()]
            elif isinstance(values, dict):
                values = [[atom.getattr(values[attr], 0.0, recursive=True)
                           for atom in ss.atoms] for ss in trajectory]
            values = np.asarray(values)
            for i, ss in enumerate(trajectory):
                ss._atoms_array[:, attridx] = values[i]

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = dict(dumpfiles=self.dumpfiles,
                         reference_timestep=self.reference_timestep,
                         reference_index=self.reference_index)
        attr_dict.update(self.formatter.todict())
        return attr_dict

DUMPIOReader = DUMPReader


class DUMPWriter:
    """Class for writing LAMMPS dump chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, structure=None,
              atoms=None, snapshot=None, trajectory=None,
              timestep=None, boxstr=None, domain=None,
              bounding_box=None, allow_triclinic_box=False,
              attr_value_map=None, **kwargs):
        """Write structure dump to file.

        Parameters
        ----------
        fname : str, optional
            Output file name.
        outpath : str, optional
            Output file path.
        fpath : str, optional
            Full path (directory path + file name) to output data file.
        atoms : `Atoms`
            An :py:class:`Atoms` instance.
        verbose : bool, optional
            verbose output

        """
        invalid_kwargs = all([obj is None for obj in
                              (atoms, structure, snapshot, trajectory)])
        if invalid_kwargs:
            error_msg = \
                ('Expected a `Trajectory` object for kwarg `trajectory`,\n'
                 '`StructureBase` sub-class for kwarg `structure`,\n'
                 'or an `Atoms` object for kwarg `atoms`.')
            raise ValueError(error_msg)

        if structure is not None and atoms is None:
            atoms = structure.atoms

        if fpath is None:
            fpath = get_fpath(fname=fname, ext='dump', outpath=outpath,
                              overwrite=True, add_fnum=False)

        dump = DUMPIO(**kwargs)
        formatter = dump.formatter
        if trajectory is None and snapshot is None:
            trajectory = dump.trajectory
            if not isinstance(atoms, MDAtoms):
                atoms = MDAtoms(atoms)
            snapshot = Snapshot(trajectory)
            if timestep is None:
                timestep = 0
            snapshot.timestep = timestep
            Natoms = snapshot.Natoms = atoms.Natoms
            snapshot.atom_selection = np.zeros(Natoms, dtype=bool)

            if boxstr is None:
                boxstr = 'pp pp pp'
            snapshot.boxstr = boxstr

            if domain is None:
                domain = snapshot.domain
                if bounding_box is not None:
                    domain.update(from_region=bounding_box,
                                  allow_triclinic_box=allow_triclinic_box,
                                  **kwargs)
                else:
                    lattice = None
                    if structure is not None and structure.lattice is not None:
                        lattice = structure.lattice
                    elif atoms.lattice is not None:
                        lattice = atoms.lattice

                    if lattice is not None:
                        domain.update(from_lattice=lattice,
                                      allow_triclinic_box=allow_triclinic_box,
                                      **kwargs)
                    else:
                        domain.update(from_array=atoms.coords,
                                      allow_triclinic_box=allow_triclinic_box,
                                      **kwargs)
            snapshot.domain = domain

            # scale_original = None
            # if all([flag is False for flag in (xflag, yflag, zflag)]):
            #     scale_original = False
            # if all([flag for flag in (xflag, yflag, zflag)]):
            #     scale_original = True

            dumpattrs = formatter.dumpattrs
            atomattrs = formatter.atomattrs
            otherattrs = formatter.otherattrs
            atoms_array = \
                np.zeros((snapshot.Natoms, len(dumpattrs)), dtype=float)
            for n, atom in enumerate(atoms):
                atoms_array[n] = [getattr(atom, attr) if attr in atomattrs
                                  else atom.getattr(attr, default=0.0,
                                                    recursive=True)
                                  if attr in otherattrs else 0.0
                                  for attr in dumpattrs]

            snapshot.atoms = atoms_array
            snapshot.formatter = formatter
            trajectory.append(snapshot)
        else:
            kwargs['atoms'] = atoms
            kwargs['snapshot'] = snapshot

        dump.write(dumpfile=fpath, **kwargs)

DUMPIOWriter = DUMPWriter


class DUMPData(DUMPReader):
    """Class for reading and writing structure data in LAMMPS dump format.

    Parameters
    ----------
    *args : one or more dump files
    **kwargs : dict mapping of keyword arguments

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def write(self, dumpfile=None, atoms=None, snapshot=None, trajectory=None,
              write_header=True, mode='w', N=None, **kwargs):
        """Write dump file.

        Parameters
        ----------
        dumpfile : {None, str}, optional
        atoms : :class:`~sknano.core.atoms.Atoms`, optional
        snapshot : :class:`~sknano.core.atoms.Snapshot`, optional
        trajectory : :class:`~sknano.core.atoms.trajectory`, optional
        write_header : :class:`~python:bool`, optional
        mode : {'w', 'a'}, optional
        N : :class:`~python:int`, optional

        """
        if trajectory is None:
            trajectory = self.trajectory

        if len(trajectory) == 0 and snapshot is not None:
            trajectory.append(snapshot)
            snapshot = None

        if trajectory.Nsnaps > 0:
            try:
                if not dumpfile:
                    if self.fpath is None:
                        error_msg = 'Invalid `dumpfile` {}'.format(dumpfile)
                        raise ValueError(error_msg)
                    else:
                        dumpfile = self.fpath
                elif self.fpath is None:
                    self.fpath = dumpfile

                trajectory.sort(key=attrgetter('timestep'))
                trajectory.cull()
                trajectory.time_selection.all()

                # self._update_attr_fmtstr_widths()

                if atoms is not None and isinstance(atoms, Atoms):
                    if not isinstance(atoms, MDAtoms):
                        atoms = MDAtoms(atoms)

                    if snapshot is None:
                        trajectory.atom_selection.update(atoms)
                    else:
                        trajectory.atom_selection.update(atoms,
                                                         ts=snapshot.timestep)

                mode += 't'

                with zopen(dumpfile, mode) as stream:
                    if snapshot is not None:
                        ss = trajectory.get_snapshot(snapshot.timestep)
                        if write_header:
                            self._write_header(stream, ss)
                            self._write_domain(stream, ss)
                        self._write_atoms(stream, ss)
                    else:
                        for ss in trajectory:
                            if not ss.selected:
                                continue
                            if write_header:
                                self._write_header(stream, ss)
                                self._write_domain(stream, ss)
                            self._write_atoms(stream, ss)

            except (OSError, TypeError, ValueError) as e:
                print(e)

    def _update_attr_fmtstr_widths(self):
        # attr_fmtstr_width = self.attr_fmtstr_width = {}
        pass

    def _write_header(self, stream, ss):
            """Write snapshot header info."""
            stream.write('ITEM: TIMESTEP\n')
            stream.write('{}\n'.format(ss.timestep))
            stream.write('ITEM: NUMBER OF ATOMS\n')
            stream.write('{}\n'.format(ss.nselected))

    def _write_domain(self, stream, ss):
            """Write snapshot bounding box info."""
            box_bounds = 'ITEM: BOX BOUNDS'
            domain = ss.domain

            if ss.boxstr:
                box_bounds = ' '.join((box_bounds, '{}'.format(ss.boxstr)))
            stream.write('{}\n'.format(box_bounds))

            for dim, tilt_factor in zip(('x', 'y', 'z'), ('xy', 'xz', 'yz')):
                lo = dim + 'lo'
                hi = dim + 'hi'
                if domain.triclinic:
                    lo = '_'.join((lo, 'bound'))
                    hi = '_'.join((hi, 'bound'))

                box_bounds = \
                    ' '.join(('{}'.format(getattr(domain, lo)),
                              '{}'.format(getattr(domain, hi))))

                if domain.triclinic:
                    box_bounds = \
                        ' '.join((box_bounds,
                                  '{}'.format(getattr(domain, tilt_factor))))
                stream.write('{}\n'.format(box_bounds))

    def _write_atoms(self, stream, ss):
        """Write snapshot atoms."""
        atoms_array = ss.get_atoms(asarray=True)[ss.atom_selection]
        formatter = ss.formatter
        attr_dtypes = formatter.attr_dtypes
        sformat = formatter.format
        stream.write('ITEM: ATOMS {}\n'.format(formatter.dumpattrs2str()))
        for atom in atoms_array:
            stream.write(sformat(atom, attr_dtypes))

DUMP = DUMPIO = DUMPData


class DUMPError(StructureDataError):
    """Exception class for :class:`DUMPData` I/O errors."""
    pass

DUMPIOError = DUMPError


class DUMPFormatter(StructureDataFormatter):
    """Class defining the structure file format for LAMMPS dump.

    Parameters
    ----------
    style : :class:`~python:str`, optional
    dumpattrs : :class:`~python:list`, optional
    dumpattrmap : class:`~python:dict`
        Python :class:`~python:dict` mapping custom dump attributes
        to :class:`~sknano.core.atoms.Atom` attributes.
    atomattrmap : :class:`~python:dict`
        :class:`~python:dict` mapping atom type to atom element

    Attributes
    ----------
    style
    dumpattrs2index

    """
    def __init__(self, style=None, dumpattrs=None, dumpattrmap=None,
                 atomattrmap=None):
        if style == 'custom' and dumpattrs is None:
            raise ValueError('Expected list of dump attrs for `dumpattrs`')

        if style is None:
            style = 'custom'
        self.style = style

        if style.startswith('atom'):
            dumpattrs = ['id', 'type', 'xs', 'ys', 'zs']

        if dumpattrs is None:
            dumpattrs = []
        # self.dumpattrs = dumpattrs

        self.dumpattrmap = dumpattrmap
        self.atomattrmap = atomattrmap

        self.dumpattrs2index = {k: i for i, k in enumerate(dumpattrs)}

        if dumpattrs:
            self._update_attrs(dumpattrs)

        self.fmtstr = "style={style!r}, " + \
            "dumpattrs={dumpattrs!r}, " + \
            "dumpattrmap={dumpattrmap!r}, " + \
            "atomattrmap={atomattrmap!r}"

    def __str__(self):
        strrep = super().__str__()
        style = self.style
        dumpattrs = self.dumpattrs
        dumpattrmap = self.dumpattrmap
        atomattrmap = self.atomattrmap
        atomattrs = self.atomattrs
        otherattrs = self.otherattrs

        items = ['style', 'dumpattrs', 'dumpattrmap', 'atomattrs',
                 'atomattrmap', 'otherattrs']
        values = [style, dumpattrs, dumpattrmap, atomattrs, atomattrmap,
                  otherattrs]
        table = self._tabulate(list(zip(items, values)))
        return '\n'.join((strrep, table))

    @property
    def attr_dtypes(self):
        """List of dump atom attributes."""
        try:
            return self._attr_dtypes
        except AttributeError:
            self._update_attr_dtypes()
            return self._attr_dtypes

    @property
    def atomattrs(self):
        """List of dump atom attributes."""
        try:
            return self._atomattrs
        except AttributeError:
            self._update_attrs()
            return self._atomattrs

    @property
    def dumpattrs(self):
        """List of dump attributes."""
        if self.dumpattrs2index:
            self._update_attrs()
            return self._dumpattrs
        else:
            return None

    @property
    def otherattrs(self):
        """List of other dump atom attributes."""
        try:
            return self._otherattrs
        except AttributeError:
            self._update_attrs()
            return self._otherattrs

    @property
    def _dumpattrslist(self):
        """List of dump attributes."""
        dumpattrs2index = self.dumpattrs2index
        return sorted(dumpattrs2index, key=dumpattrs2index.__getitem__)

    def _update_attrs(self, dumpattrs=None):
        """Update :class:`~python:dict` mapping of unknown dump \
            attributes->:attr:`~DUMPFormatter.dumpattrs` list index."""
        if dumpattrs is None:
            dumpattrs = self._dumpattrslist
        dumpattrmap = self.dumpattrmap
        if dumpattrmap is not None:
            self.remapdumpattrs(dumpattrs, dumpattrmap)
        self._dumpattrs = dumpattrs[:]
        atomattrs = dumpattrs[:]
        otherattrmap = {attr: dumpattrs.index(attr) for attr in
                        set(dumpattrs) - set(dir(Atom()))}
        otherattrs = sorted(otherattrmap, key=otherattrmap.__getitem__)
        [atomattrs.remove(attr) for attr in otherattrs]
        self._atomattrs = atomattrs
        self._otherattrs = otherattrs

    def _update_attr_dtypes(self):
        self._attr_dtypes = [DUMPATTR_DTYPES[attr] if attr in DUMPATTR_DTYPES
                             else float for attr in self.dumpattrs]

    def atomattrs2str(self):
        """Return a space-separated string of parsed atom attributes."""
        return ' '.join(self.atomattrs)

    def dumpattrs2str(self):
        """Return a space-separated string of the dumped atom attributes."""
        return ' '.join(self.dumpattrs)

    def remapdumpattrs(self, dumpattrs, dumpattrmap):
        """Rename attributes in the :attr:`DUMPFormatter.dumpattrs` list.

        Parameters
        ----------
        dumpattrmap : :class:`~python:dict`
            :class:`~python:dict` mapping atom attribute name to new
            attribute name.

        """
        for k, v in dumpattrmap.items():
            try:
                dumpattrs[dumpattrs.index(k)] = v
            except ValueError:
                pass

    def format(self, attrs, attr_dtypes=None):
        """Return :class:`~python:str` of dump attributes formatted for an \
            output stream.

        Parameters
        ----------
        attrs : :class:`~python:list`
            List of dump attribute values.

        """
        if attr_dtypes is None:
            attr_dtypes = self.attr_dtypes
        line = ' '.join(('{}'.format(dtype(attr)) for attr, dtype in
                        zip(attrs, attr_dtypes)))
        line += '\n'
        return line

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(style=self.style,
                              dumpattrs=self.dumpattrs,
                              dumpattrmap=self.dumpattrmap,
                              atomattrmap=self.atomattrmap))
        return attr_dict

DUMPIOFormatter = DUMPFormatter

DUMPATTR_DTYPES = \
    {'id': int, 'type': int, 'mol': int, 'q': float, 'mass': float,
     'x': float, 'y': float, 'z': float,
     'ix': int, 'iy': int, 'iz': int,
     'vx': float, 'vy': float, 'vz': float,
     'fx': float, 'fy': float, 'fz': float,
     'lx': float, 'ly': float, 'lz': float,
     'wx': float, 'wy': float, 'wz': float,
     'shapex': float, 'shapey': float, 'shapez': float,
     'quatw': float, 'quati': float, 'quatj': float, 'quatk': float,
     'atom1': int, 'atom2': int, 'atom3': int, 'atom4': int}
