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
import warnings

import numpy as np

from monty.io import zopen
from sknano.core import get_fpath, flatten
from sknano.core.atoms import Trajectory, Snapshot, Atoms, MDAtoms, \
    MDAtom as Atom
# from sknano.core.crystallography import Crystal3DLattice
from sknano.core.geometric_regions import Cuboid

from ._base import StructureIO, StructureIOError, StructureFormatSpec, \
    default_comment_line

__all__ = ['DUMPData', 'DUMPReader', 'DUMPWriter', 'DUMPIOError',
           'DUMPFormatSpec']

attr_dtypes = {'id': int, 'type': int, 'mol': int, 'q': float, 'mass': float,
               'x': float, 'y': float, 'z': float,
               'ix': int, 'iy': int, 'iz': int,
               'vx': float, 'vy': float, 'vz': float,
               'fx': float, 'fy': float, 'fz': float,
               'lx': float, 'ly': float, 'lz': float,
               'wx': float, 'wy': float, 'wz': float,
               'shapex': float, 'shapey': float, 'shapez': float,
               'quatw': float, 'quati': float, 'quatj': float, 'quatk': float,
               'atom1': int, 'atom2': int, 'atom3': int, 'atom4': int}


# class LAMMPSBOX(Crystal3DLattice):
#     """LAMMPS 3D simulation box.

#     Parameters
#     ----------
#     lx, ly, lz : float
#     xy, xz, yz : float

#     """

#     def __init__(self, lx=None, ly=None, lz=None,
#                  xy=None, xz=None, yz=None, cell_matrix=None,
#                  orientation_matrix=None, offset=None):

#         a = lx
#         b = np.sqrt(ly ** 2 + xy ** 2)
#         c = np.sqrt(lz ** 2 + xz ** 2 + yz ** 2)
#         alpha = np.degrees(np.arccos((xy * xz + ly * yz) / (b * c)))
#         beta = np.degrees(np.arccos(xz / c))
#         gamma = np.degrees(np.arccos(xy / b))
#         super().__init__(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,
#                          offset=offset)


class DUMPReader(StructureIO):
    """Class for reading `LAMMPS dump` file format.

    Parameters
    ----------
    *args : :class:`~python:list`
        :class:`~python:list` of one or more LAMMPS dump files.
    dumpattrmap : class:`~python:dict`
        Python :class:`~python:dict` mapping custom dump attributes
        to :class:`~sknano.core.atoms.Atom` attributes.
    elementmap : :class:`~python:dict`
        :class:`~python:dict` mapping atom type to atom element
    autoread : :class:`~python:bool`, optional
        Automatically read dump files. Default is `True`

        .. versionadded:: 0.3.22

    reference_timestep : :class:`~python:int`, optional
        The `timestep` corresponding to the
        :attr:`~sknano.core.atoms.Snapshot.timestep` of the
        :class:`~sknano.core.atoms.Trajectory`\ s
        :attr:`~sknano.core.atoms.Trajectory.reference_snapshot`.
        Default is `None`.

        .. note::
           **Overrides** `reference_index` if not `None`

        .. versionadded:: 0.3.22

    reference_index : :class:`~python:int`, optional
        The :class:`~python:list` index corresponding to the
        :attr:`~sknano.core.atoms.Snapshot.timestep` of the
        :class:`~sknano.core.atoms.Trajectory`\ s
        :attr:`~sknano.core.atoms.Trajectory.reference_snapshot`.
        Default is `None`.

        .. versionadded:: 0.3.22

    Attributes
    ----------
    trajectory
    dumpattrs
    atomattrs

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
    The `elementmap` parameter takes a :class:`~python:dict` mapping
    of the LAMMPS dump `type` value to an element symbol::

    >>> dumps = DUMPReader('dump.*', dumpattrmap=dumpattrmap,
    ...                    elementmap={1: 'C'})
    >>> atoms = dumps[0].atoms
    >>> print(atoms[0])

    """
    def __init__(self, *args, autoread=True, dumpattrmap=None,
                 elementmap=None, reference_timestep=None,
                 reference_index=None, **kwargs):

        if 'attrmap' in kwargs:
            msg = ("The {!s} keyword argument `attrmap` was deprecated in "
                   "in version 0.4.0.\nUse `dumpattrmap` instead."
                   .format(self.__class__.__name__))
            warnings.warn(msg, DeprecationWarning, stacklevel=2)
            dumpattrmap = kwargs.pop('attrmap')

        super().__init__(**kwargs)

        self.dumpattrmap = dumpattrmap
        self.elementmap = elementmap
        self.trajectory = Trajectory()
        self.dumpattrs = {}
        self.dumpfiles = tuple(flatten([[glob(f) for f in arg.split()]
                                        for arg in args]))
        self._reference_timestep = reference_timestep
        self._reference_index = reference_index

        if len(self.dumpfiles) == 0:
            raise ValueError('No dump file specified.')

        self.fmtstr = "{dumpfiles!r}, autoread=True, " + \
            "dumpattrmap={dumpattrmap!r}, " + \
            "elementmap={elementmap!r}, " + \
            "reference_timestep={reference_timestep!r}, " + \
            "reference_index={reference_index!r}"

        if autoread:
            self.read()

    def __getattr__(self, name):
        try:
            return getattr(self.trajectory, name)
        except AttributeError:
            return super().__getattr__(name)

    def __iter__(self):
        return iter(self.trajectory)

    def __getitem__(self, index):
        return self.trajectory[index]

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
        for dumpfile in self.dumpfiles:
            with zopen(dumpfile) as f:
                snapshot = self.read_snapshot(f)
                while snapshot is not None:
                    self.trajectory.append(snapshot)
                    print(snapshot.timestep, end=' ')
                    sys.stdout.flush()
                    snapshot = self.read_snapshot(f)
        print()

        self.trajectory.sort(key=attrgetter('timestep'))
        self.trajectory.cull()

        print("read {:d} snapshots".format(self.Nsnaps))

        self.trajectory.time_selection.all()
        self._update_reference_snapshot()

        if self.dumpattrs:
            print('Dumped Atom attributes: {}'.format(self.dumpattrs2str()))
        else:
            print('No dump column assignments')

        if 'x' not in self.dumpattrs or 'y' not in self.dumpattrs or \
                'z' not in self.dumpattrs:
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
            f.readline()
            snapshot.timestep = int(f.readline().strip().split()[0])
            f.readline()
            snapshot.Natoms = int(f.readline().strip())
            snapshot.atom_selection = np.zeros(snapshot.Natoms, dtype=bool)

            item = f.readline().strip()
            try:
                snapshot.boxstr = item.split('BOUNDS')[1].strip()
            except IndexError:
                snapshot.boxstr = ''

            snapshot.triclinic = False
            snapshot.bounding_box = Cuboid()

            if 'xy' in snapshot.boxstr:
                snapshot.triclinic = True

            for dim, tilt_factor in zip(('x', 'y', 'z'), ('xy', 'xz', 'yz')):
                bounds = f.readline().strip().split()
                setattr(snapshot, dim + 'lo', float(bounds[0]))
                setattr(snapshot, dim + 'hi', float(bounds[1]))
                setattr(snapshot.bounding_box, dim + 'min',
                        getattr(snapshot, dim + 'lo'))
                setattr(snapshot.bounding_box, dim + 'max',
                        getattr(snapshot, dim + 'hi'))
                try:
                    setattr(snapshot, tilt_factor, float(bounds[2]))
                except IndexError:
                    setattr(snapshot, tilt_factor, 0.0)

            if not self.dumpattrs:
                xflag = yflag = zflag = None
                attrs = f.readline().strip().split()[2:]
                for i, attr in enumerate(attrs):
                    if attr in ('x', 'xu', 'xs', 'xsu'):
                        self.dumpattrs['x'] = i
                        if attr in ('x', 'xu'):
                            xflag = False
                        else:
                            xflag = True
                    elif attr in ('y', 'yu', 'ys', 'ysu'):
                        self.dumpattrs['y'] = i
                        if attr in ('y', 'yu'):
                            yflag = False
                        else:
                            yflag = True
                    elif attr in ('z', 'zu', 'zs', 'zsu'):
                        self.dumpattrs['z'] = i
                        if attr in ('z', 'zu'):
                            zflag = False
                        else:
                            zflag = True
                    else:
                        self.dumpattrs[attr] = i

                self.scale_original = None
                if all([flag is False for flag in (xflag, yflag, zflag)]):
                    self.scale_original = False
                if all([flag for flag in (xflag, yflag, zflag)]):
                    self.scale_original = True

                self.atomattrs = \
                    sorted(self.dumpattrs, key=self.dumpattrs.__getitem__)

                self.attr_dtypes = [attr_dtypes[attr] if attr in attr_dtypes
                                    else float for attr in self.atomattrs]

                if self.dumpattrmap is not None:
                    self.remap_atomattr_names(self.dumpattrmap)
                    self.dumpattrmap = None

                self.unknown_attrs = \
                    {attr: self.atomattrs.index(attr) for
                     attr in set(self.atomattrs) - set(dir(Atom()))}

                [self.atomattrs.remove(attr) for attr in self.unknown_attrs]
            else:
                f.readline()

            snapshot.atomattrs = self.atomattrs
            snapshot.attr_dtypes = self.attr_dtypes
            snapshot.elementmap = self.elementmap

            atoms = \
                np.zeros((snapshot.Natoms, len(self.atomattrs)), dtype=float)
            for n in range(snapshot.Natoms):
                line = [float(value) for col, value in
                        enumerate(f.readline().strip().split())
                        if col not in self.unknown_attrs.values()]
                atoms[n] = line

            snapshot.atoms = atoms
            return snapshot

        except IndexError:
            return None

    def remap_atomattr_names(self, dumpattrmap):
        """Rename attributes in the :attr:`DUMPReader.atomattrs` list.

        Parameters
        ----------
        dumpattrmap : :class:`~python:dict`
            :class:`~python:dict` mapping atom attribute name to new
            attribute name.

        """
        for k, v in dumpattrmap.items():
            try:
                self.atomattrs[self.atomattrs.index(k)] = v
            except ValueError:
                pass

    def scale(self):
        """Scale cartesian coordinates to fractional coordinates."""
        xi = self.dumpattrs['x']
        yi = self.dumpattrs['y']
        zi = self.dumpattrs['z']
        for snapshot in self.trajectory:
            atoms = snapshot.get_atoms(asarray=True)
            if atoms is not None:
                xlo = snapshot.xlo
                xhi = snapshot.xhi
                ylo = snapshot.ylo
                yhi = snapshot.yhi
                zlo = snapshot.zlo
                zhi = snapshot.zhi
                lx = xhi - xlo
                ly = yhi - ylo
                lz = zhi - zlo
                xy = snapshot.xy
                xz = snapshot.xz
                yz = snapshot.yz

                if np.allclose([xy, xz, yz], np.zeros(3)):
                    atoms[:, xi] = (atoms[:, xi] - snapshot.xlo) / lx
                    atoms[:, yi] = (atoms[:, yi] - snapshot.ylo) / ly
                    atoms[:, zi] = (atoms[:, zi] - snapshot.zlo) / lz
                else:
                    xlo = xlo - min((0.0, xy, xz, xy + xz))
                    xhi = xhi - max((0.0, xy, xz, xy + xz))
                    lx = xhi - xlo

                    ylo = ylo - min((0.0, yz))
                    yhi = yhi - max((0.0, yz))
                    ly = yhi - ylo

                    atoms[:, xi] = (atoms[:, xi] - snapshot.xlo) / lx + \
                        (atoms[:, yi] - snapshot.ylo) * xy / (lx * ly) + \
                        (atoms[:, zi] - snapshot.zlo) * (yz * xy - ly * xz) / \
                        (lx * ly * lz)
                    atoms[:, yi] = (atoms[:, yi] - snapshot.ylo) / ly + \
                        (atoms[:, zi] - snapshot.zlo) * yz / (ly * lz)
                    atoms[:, zi] = (atoms[:, zi] - snapshot.zlo) / lz

    def unscale(self):
        """Unscale fractional coordinates to cartesian coordinates."""
        xi = self.dumpattrs['x']
        yi = self.dumpattrs['y']
        zi = self.dumpattrs['z']
        for snapshot in self.trajectory:
            atoms = snapshot.get_atoms(asarray=True)
            if atoms is not None:
                xlo = snapshot.xlo
                xhi = snapshot.xhi
                ylo = snapshot.ylo
                yhi = snapshot.yhi
                zlo = snapshot.zlo
                zhi = snapshot.zhi
                lx = xhi - xlo
                ly = yhi - ylo
                lz = zhi - zlo
                xy = snapshot.xy
                xz = snapshot.xz
                yz = snapshot.yz
                if np.allclose([xy, xz, yz], np.zeros(3)):
                    atoms[:, xi] = snapshot.xlo + atoms[:, xi] * lx
                    atoms[:, yi] = snapshot.ylo + atoms[:, yi] * ly
                    atoms[:, zi] = snapshot.zlo + atoms[:, zi] * lz
                else:
                    xlo = xlo - min((0.0, xy, xz, xy + xz))
                    xhi = xhi - max((0.0, xy, xz, xy + xz))
                    lx = xhi - xlo

                    ylo = ylo - min((0.0, yz))
                    yhi = yhi - max((0.0, yz))
                    ly = yhi - ylo

                    atoms[:, xi] = snapshot.xlo + atoms[:, xi] * lx + \
                        atoms[:, yi] * xy + atoms[:, zi] * xz
                    atoms[:, yi] = snapshot.ylo + atoms[:, yi] * ly + \
                        atoms[:, zi] * yz
                    atoms[:, zi] = snapshot.zlo + atoms[:, zi] * lz

    def atomattrs2str(self):
        """Return a space-separated string of parsed atom attributes."""
        return ' '.join(self.atomattrs)

    def dumpattrs2str(self):
        """Return a space-separated string of the dumped atom attributes."""
        return ' '.join(sorted(self.dumpattrs, key=self.dumpattrs.__getitem__))

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(dumpfiles=self.dumpfiles, dumpattrmap=self.dumpattrmap,
                    elementmap=self.elementmap,
                    reference_timestep=self.reference_timestep,
                    reference_index=self.reference_index)


class DUMPWriter:
    """Class for writing LAMMPS dump chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, atoms=None,
              comment_line=None, verbose=False, **kwargs):
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
        boxbounds : dict, optional
            If `None`, determined automatically from atom coordinates.
        comment_line : str, optional
            A string written to the first line of `dump` file. If `None`,
            then it is set to the full path of the output `dump` file.
        assert_unique_ids : bool, optional
            Check that each Atom in Atoms has a unique ID. If the check
            fails, then assign a unique ID to each Atom.
            If `assert_unique_ids` is True, but the ID's are not
            unique, LAMMPS will not be able to read the dump file.
        verbose : bool, optional
            verbose output

        """
        if fpath is None:
            fpath = get_fpath(fname=fname, ext='dump', outpath=outpath,
                              overwrite=True, add_fnum=False)

        if comment_line is None:
            comment_line = default_comment_line

        atoms.rezero_coords()


class DUMPData(DUMPReader):
    """Class for reading and writing structure data in LAMMPS dump format.

    Parameters
    ----------
    *args : one or more dump files

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def write(self, dumpfile=None, snapshot=None, atoms=None, dumpattrs=None):
        """Write dump file.

        Parameters
        ----------
        dumpfile : {None, str}, optional

        """
        if self.trajectory.Nsnaps > 0:
            try:
                if not dumpfile:
                    if self.fpath is None:
                        error_msg = 'Invalid `dumpfile` {}'.format(dumpfile)
                        raise ValueError(error_msg)
                    else:
                        dumpfile = self.fpath
                # DUMPWriter.write(fname=dumpfile, atoms=self.atoms,
                #                  comment_line=self.comment_line)

                # self._update_attr_fmtstr_widths()

                if atoms is not None and isinstance(atoms, Atoms):
                    if not isinstance(atoms, MDAtoms):
                        atoms = MDAtoms(atoms)

                    if snapshot is None:
                        self.trajectory.atom_selection.update(atoms)
                    else:
                        self.trajectory.atom_selection.update(
                            atoms, ts=snapshot.timestep)

                with zopen(dumpfile, 'wt') as fh:
                    if snapshot is not None:
                        ss = self.trajectory.get_snapshot(snapshot.timestep)
                        self._write_header(fh, ss)
                        self._write_bounding_box(fh, ss)
                        self._write_atoms(fh, ss)
                    else:
                        for ss in self.trajectory:
                            if not ss.selected:
                                continue
                            self._write_header(fh, ss)
                            self._write_bounding_box(fh, ss)
                            self._write_atoms(fh, ss)

            except (OSError, TypeError, ValueError) as e:
                print(e)

    def _update_attr_fmtstr_widths(self):
        attr_fmtstr_width = self.attr_fmtstr_width = {}
        # for attr in self.atomattrs:
        #     attr_fmtstr_width[attr] = len(str(self))

        # attr_fmtstr_width['id'] = len(str(self.atoms.Natoms)) + 1
        # attr_fmtstr_width['type'] = len(str(self.atoms.Ntypes)) + 1
        # attr_fmtstr_width['mol'] = len(str(np.max(self.atoms.mols))) + 1
        # attr_fmtstr_width['q'] = \
        #     len('{:f}'.format(np.max(self.atoms.charges))) + 2
        # attr_fmtstr_width['mass'] = \
        #     len('{:f}'.format(np.max(self.atoms.masses))) + 4
        # attr_fmtstr_width['x'] = len('{:f}'.format(np.max(self.atoms.x))) + 2
        # attr_fmtstr_width['y'] = len('{:f}'.format(np.max(self.atoms.y))) + 2
        # attr_fmtstr_width['z'] = len('{:f}'.format(np.max(self.atoms.z))) + 2
        # attr_fmtstr_width['ix'] = len(str(np.max(self.atoms.ix))) + 2
        # attr_fmtstr_width['iy'] = len(str(np.max(self.atoms.iy))) + 2
        # attr_fmtstr_width['iz'] = len(str(np.max(self.atoms.iz))) + 2
        # attr_fmtstr_width['vx'] = len('{:f}'.format(np.max(self.atoms.vx))) + 2
        # attr_fmtstr_width['vy'] = len('{:f}'.format(np.max(self.atoms.vy))) + 2
        # attr_fmtstr_width['vz'] = len('{:f}'.format(np.max(self.atoms.vz))) + 2

    def _write_header(self, fh, ss):
            """Write snapshot header info."""
            fh.write('ITEM: TIMESTEP\n')
            fh.write('{}\n'.format(ss.timestep))
            fh.write('ITEM: NUMBER OF ATOMS\n')
            fh.write('{}\n'.format(ss.nselected))

    def _write_bounding_box(self, fh, ss):
            """Write snapshot bounding box info."""
            box_bounds = 'ITEM: BOX BOUNDS'
            if ss.boxstr:
                box_bounds = ' '.join((box_bounds, '{}'.format(ss.boxstr)))
            fh.write('{}\n'.format(box_bounds))
            for dim, tilt_factor in zip(('x', 'y', 'z'), ('xy', 'xz', 'yz')):
                box_bounds = ' '.join(('{}'.format(getattr(ss, dim + 'lo')),
                                       '{}'.format(getattr(ss, dim + 'hi'))))
                if ss.triclinic:
                    box_bounds = \
                        ' '.join((box_bounds,
                                  '{}'.format(getattr(ss, tilt_factor))))
                fh.write('{}\n'.format(box_bounds))

    def _write_atoms(self, fh, ss):
        """Write snapshot atoms."""
        fh.write('ITEM: ATOMS {}\n'.format(self.atomattrs2str()))
        aselection = ss.get_atoms(asarray=True)[ss.atom_selection]
        attr_dtypes = ss.attr_dtypes
        for atom in aselection:
            line = ' '.join(('{}'.format(dtype(attr)) for dtype, attr in
                            zip(attr_dtypes, atom)))
            line += '\n'
            fh.write(line)


class DUMPIOError(StructureIOError):
    """Exception class for :class:`DUMPData` I/O errors."""
    pass


class DUMPFormatSpec(StructureFormatSpec):
    """Class defining the structure file format for LAMMPS dump."""
    pass
