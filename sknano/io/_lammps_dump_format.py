# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS dump format (:mod:`sknano.io._lammps_dump_format`)
====================================================================

.. currentmodule:: sknano.io._lammps_dump_format

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import glob
# import re
import sys
from operator import attrgetter

import numpy as np

from monty.io import zopen
from sknano.core import get_fpath
from sknano.core.atoms import Trajectory, Snapshot, MDAtom as Atom
from sknano.core.crystallography import Crystal3DLattice
from sknano.core.geometric_regions import Cuboid, Parallelepiped

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


class DUMPReader(StructureIO):
    """Class for reading `LAMMPS dump` file format.

    Parameters
    ----------
    fpath : :class:`~python:str`
        LAMMPS dump file path

    """
    def __init__(self, *args, attrmap=None, **kwargs):
        super().__init__(**kwargs)

        self.attrmap = attrmap
        self.trajectory = Trajectory()
        self.dumpattrs = {}
        self.dumpfiles = []
        for fpath in args:
            self.dumpfiles.append(fpath)

        if len(self.dumpfiles) == 0:
            raise ValueError('No dump file specified.')

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
        self.trajectory.t0_snapshot = self.trajectory[0]

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

                if self.attrmap is not None:
                    self.remap_atomattr_names(self.attrmap)
                    self.attrmap = None

                self.unknown_attrs = \
                    {attr: self.atomattrs.index(attr) for
                     attr in set(self.atomattrs) - set(dir(Atom()))}

                [self.atomattrs.remove(attr) for attr in self.unknown_attrs]
            else:
                f.readline()

            snapshot.atomattrs = self.atomattrs
            snapshot.attr_dtypes = self.attr_dtypes

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

    def remap_atomattr_names(self, attrmap):
        """Rename attributes in the :attr:`DUMPReader.atomattrs` list.

        Parameters
        ----------
        attrmap : :class:`~python:dict`
            :class:`~python:dict` mapping atom attribute name to new
            attribute name.

        """
        for k, v in attrmap.items():
            try:
                self.atomattrs[self.atomattrs.index(k)] = v
            except ValueError:
                pass

    def scale(self):
        x = self.dumpattrs['x']
        y = self.dumpattrs['y']
        z = self.dumpattrs['z']
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
                    atoms[:, x] = (atoms[:, x] - snapshot.xlo) / lx
                    atoms[:, y] = (atoms[:, y] - snapshot.ylo) / ly
                    atoms[:, z] = (atoms[:, z] - snapshot.zlo) / lz
                else:
                    xlo = xlo - min((0.0, xy, xz, xy + xz))
                    xhi = xhi - max((0.0, xy, xz, xy + xz))
                    lx = xhi - xlo

                    ylo = ylo - min((0.0, yz))
                    yhi = yhi - max((0.0, yz))
                    ly = yhi - ylo

                    atoms[:, x] = (atoms[:, x] - snapshot.xlo) / lx + \
                        (atoms[:, y] - snapshot.ylo) * xy / (lx * ly) + \
                        (atoms[:, z] - snapshot.zlo) * (yz * xy - ly * xz) / \
                        (lx * ly * lz)
                    atoms[:, y] = (atoms[:, y] - snapshot.ylo) / ly + \
                        (atoms[:, z] - snapshot.zlo) * yz / (ly * lz)
                    atoms[:, z] = (atoms[:, z] - snapshot.zlo) / lz

    def unscale(self):
        x = self.dumpattrs['x']
        y = self.dumpattrs['y']
        z = self.dumpattrs['z']
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
                    atoms[:, x] = snapshot.xlo + atoms[:, x] * lx
                    atoms[:, y] = snapshot.ylo + atoms[:, y] * ly
                    atoms[:, z] = snapshot.zlo + atoms[:, z] * lz
                else:
                    xlo = xlo - min((0.0, xy, xz, xy + xz))
                    xhi = xhi - max((0.0, xy, xz, xy + xz))
                    lx = xhi - xlo

                    ylo = ylo - min((0.0, yz))
                    yhi = yhi - max((0.0, yz))
                    ly = yhi - ylo

                    atoms[:, x] = snapshot.xlo + \
                        atoms[:, x] * lx + atoms[:, y] * xy + atoms[:, z] * xz
                    atoms[:, y] = snapshot.ylo + \
                        atoms[:, y] * ly + atoms[:, z] * yz
                    atoms[:, z] = snapshot.zlo + atoms[:, z] * lz

    def dumpattrs2str(self):
        return ' '.join(sorted(self.dumpattrs, key=self.dumpattrs.__getitem__))


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
    def __init__(self, *args):
        super().__init__(*args)

    def write(self, dumpfile=None):
        """Write dump file.

        Parameters
        ----------
        dumpfile : {None, str}, optional

        """
        try:
            if (dumpfile is None or dumpfile == '') and \
                    (self.fpath is None or self.fpath == ''):
                error_msg = '`dumpfile` must be a string at least 1 ' + \
                    'character long.'
                if dumpfile is None:
                    raise TypeError(error_msg)
                else:
                    raise ValueError(error_msg)
            else:
                dumpfile = self.fpath
            DUMPWriter.write(fname=dumpfile, atoms=self.atoms,
                             comment_line=self.comment_line)
        except (TypeError, ValueError) as e:
            print(e)


class DUMPIOError(StructureIOError):
    pass


class DUMPFormatSpec(StructureFormatSpec):
    """Class defining the structure file format for LAMMPS dump."""
    pass
