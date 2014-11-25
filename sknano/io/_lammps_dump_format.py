# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS dump format (:mod:`sknano.io._lammps_dump_format`)
====================================================================

.. currentmodule:: sknano.io._lammps_dump_format

"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip
__docformat__ = 'restructuredtext en'

#import glob
#import re
import sys
from operator import attrgetter

import numpy as np

from sknano.core import get_fpath
from sknano.core.atoms import Trajectory, Snapshot

from ._base import StructureIO, StructureIOError, StructureFormatSpec, \
    default_comment_line

__all__ = ['DUMPData', 'DUMPReader', 'DUMPWriter', 'DUMPIOError',
           'DUMPFormatSpec']

attr_dtypes = {'atomID': int, 'atomtype': int, 'bondID': int, 'bondtype': int,
               'moleculeID': int, 'q': float, 'ervel': float,
               'm': float, 'mass': float,
               'x': float, 'y': float, 'z': float,
               'nx': int, 'ny': int, 'nz': int,
               'vx': float, 'vy': float, 'vz': float,
               'lx': float, 'ly': float, 'lz': float,
               'wx': float, 'wy': float, 'wz': float,
               'fx': float, 'fy': float, 'fz': float,
               'atom_id': int, 'molecule_id': int, 'bond_id': int,
               'atom1': int, 'atom2': int, 'atom3': int, 'atom4': int,
               'dihedral_id': int, 'dihedraltype': int,
               'shapex': float, 'shapey': float, 'shapez': float,
               'quatw': float, 'quati': float, 'quatj': float, 'quatk': float}


class DUMPReader(StructureIO):
    """Class for reading `LAMMPS dump` file format.

    Parameters
    ----------
    fpath : str
        LAMMPS dump file path

    """
    def __init__(self, *args, **kwargs):
        super(DUMPReader, self).__init__(**kwargs)

        self.trajectory = Trajectory()
        self.dumpattrs = {}
        self.dumpfiles = []
        for fpath in args:
            self.dumpfiles.append(fpath)

        if len(self.dumpfiles) == 0:
            raise ValueError('No dump file specified.')

        self.read_all()

    @property
    def Nsnaps(self):
        return self.trajectory.Nsnaps

    def read_all(self):
        """Read all snapshots from each dump file."""
        for dumpfile in self.dumpfiles:
            with open(dumpfile) as f:
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

        self.trajectory.tselect.all()
        if self.dumpattrs:
            print('Dumped Atom attributes: {}'.format(self.dumpattrs2str()))
        else:
            print('No dump column assignments')

        if 'x' not in self.dumpattrs or 'y' not in self.dumpattrs or \
                'z' not in self.dumpattrs:
            print('dump scaling status unknown')
        elif self.Nsnaps > 0:
            if self.scale_original == 1:
                self.unscale()
            elif self.scale_original == 0:
                print('dump is already unscaled')
            else:
                print('dump scaling status unknown')

    def read_snapshot(self, f):
        try:
            snapshot = Snapshot()
            f.readline()
            snapshot.timestep = int(f.readline().strip().split()[0])
            f.readline()
            snapshot.Natoms = int(f.readline().strip())
            snapshot.aselect = np.zeros(snapshot.Natoms)

            item = f.readline().strip()
            try:
                snapshot.boxstr = item.split('BOUNDS')[1].strip()
            except IndexError:
                snapshot.boxstr = ''

            snapshot.triclinic = False
            if 'xy' in snapshot.boxstr:
                snapshot.triclinic = True

            for axis, sf in zip(('x', 'y', 'z'), ('xy', 'xz', 'yz')):
                bounds = f.readline().strip().split()
                setattr(snapshot, axis + 'lo', float(bounds[0]))
                setattr(snapshot, axis + 'hi', float(bounds[1]))
                try:
                    setattr(snapshot, sf, float(bounds[-1]))
                except IndexError:
                    setattr(snapshot, sf, 0.0)

            if self.dumpattrs:
                f.readline()
            else:
                self.scale_original = -1
                xflag = yflag = zflag = -1
                attrs = f.readline().strip().split()[2:]
                for i, attr in enumerate(attrs):
                    if attr in ('x', 'xu', 'xs', 'xsu'):
                        self.dumpattrs['x'] = i
                        if attr in ('x', 'xu'):
                            xflag = 0
                        else:
                            xflag = 1
                    elif attr in ('y', 'yu', 'ys', 'ysu'):
                        self.dumpattrs['y'] = i
                        if attr in ('y', 'yu'):
                            yflag = 0
                        else:
                            yflag = 1
                    elif attr in ('z', 'zu', 'zs', 'zsu'):
                        self.dumpattrs['z'] = i
                        if attr in ('z', 'zu'):
                            zflag = 0
                        else:
                            zflag = 1
                    else:
                        self.dumpattrs[attr] = i
                if xflag == yflag == zflag == 0:
                    self.scale_original = 0
                if xflag == yflag == zflag == 1:
                    self.scale_original = 1

                self.atomattrs = \
                    sorted(self.dumpattrs, key=self.dumpattrs.__getitem__)
                if 'id' in self.atomattrs:
                    self.atomattrs[self.atomattrs.index('id')] = 'atomID'

                if 'type' in self.atomattrs:
                    self.atomattrs[self.atomattrs.index('type')] = 'atomtype'

                self.attr_dtypes = [attr_dtypes[attr] if attr in attr_dtypes
                                    else float for attr in self.atomattrs]

            snapshot.atomattrs = self.atomattrs
            snapshot.attr_dtypes = self.attr_dtypes

            atoms = \
                np.zeros((snapshot.Natoms, len(self.atomattrs)), dtype=float)
            for n in range(snapshot.Natoms):
                line = [float(attr) for attr in f.readline().strip().split()]
                atoms[n] = line

            snapshot.atoms = atoms
            return snapshot

        except IndexError:
            return None

    def next(self):
        """Read next snapshot from list of dump files."""
        #if not self.increment:
        #    raise DUMPIOError('cannot read dump incrementally')

        #while True:
        #    with open(self.dumpfiles[self.nextfile]) as f:
        #        f.seek
        pass

    def scale(self, ts=None):
        if ts is None:
            x = self.dumpattrs['x']
            y = self.dumpattrs['y']
            z = self.dumpattrs['z']
            for snapshot in self.trajectory:
                self.scale_one(snapshot, x, y, z)

    def scale_one(self, snap, x, y, z):
        pass

    def unscale(self, ts=None):
        if ts is None:
            x = self.dumpattrs['x']
            y = self.dumpattrs['y']
            z = self.dumpattrs['z']
            for snapshot in self.trajectory:
                self.unscale_one(snapshot, x, y, z)

    def unscale_one(self, snapshot, x, y, z):
        if snapshot.xy == snapshot.xz == snapshot.yz == 0.0:
            xprd = snapshot.xhi - snapshot.xlo
            yprd = snapshot.yhi - snapshot.ylo
            zprd = snapshot.zhi - snapshot.zlo
            atoms = snapshot.get_atoms(asarray=True)
            if atoms is not None:
                atoms[:, x] = snapshot.xlo + atoms[:, x] * xprd
                atoms[:, y] = snapshot.ylo + atoms[:, y] * yprd
                atoms[:, z] = snapshot.zlo + atoms[:, z] * zprd
        else:
            xlo_bound = snapshot.xlo
            xhi_bound = snapshot.xhi
            ylo_bound = snapshot.ylo
            yhi_bound = snapshot.yhi
            zlo_bound = snapshot.zlo
            zhi_bound = snapshot.zhi
            xy = snapshot.xy
            xz = snapshot.xz
            yz = snapshot.yz
            xlo = xlo_bound - min((0.0, xy, xz, xy + xz))
            xhi = xhi_bound - max((0.0, xy, xz, xy + xz))
            ylo = ylo_bound - min((0.0, yz))
            yhi = yhi_bound - max((0.0, yz))
            zlo = zlo_bound
            zhi = zhi_bound
            h0 = xhi - xlo
            h1 = yhi - ylo
            h2 = zhi - zlo
            h3 = yz
            h4 = xz
            h5 = xy
            atoms = snapshot.get_atoms(asarray=True)
            if atoms is not None:
                atoms[:, x] = snapshot.xlo + \
                    atoms[:, x] * h0 + atoms[:, y] * h5 + atoms[:, z] * h4
                atoms[:, y] = snapshot.ylo + \
                    atoms[:, y] * h1 + atoms[:, z] * h3
                atoms[:, z] = snapshot.zlo + atoms[:, z] * h2

    def wrap(self):
        pass

    def unwrap(self):
        pass

    def owrap(self):
        pass

    def atom(self, n, fields=None):
        pass

    def dumpattrs2str(self):
        #return self.dumpattrs
        return ' '.join(sorted(self.dumpattrs, key=self.dumpattrs.__getitem__))

    def sort(self, *list):
        pass

    def sort_one(self, snapshot, id):
        pass

    def scatter(self, root):
        pass

    def minmax(self, colname):
        pass

    def set(self, eq):
        pass

    def setv(self, colname, vec):
        pass

    def clone(self, nstep, col):
        pass

    def spread(self, old, n, new):
        pass

    def vecs(self, n, *list):
        pass

    def newcolumn(self, str):
        pass

    def delete(self):
        pass

    def extra(self):
        pass

    def maxbox(self):
        pass

    def maxtype(self):
        pass


class DUMPWriter(object):
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
        super(DUMPData, self).__init__(*args)

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
