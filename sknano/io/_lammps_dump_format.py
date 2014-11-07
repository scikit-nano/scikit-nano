# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS dump format (:mod:`sknano.io._lammps_dump_format`)
====================================================================

.. currentmodule:: sknano.io._lammps_dump_format

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import glob
#import re
import sys
from operator import attrgetter

import numpy as np

from sknano.core import get_fpath
#from sknano.core.atoms import Trajectory

from ._base import Atom, Atoms, \
    StructureIO, StructureIOError, StructureFormatSpec, \
    default_comment_line

__all__ = ['DUMPData', 'DUMPReader', 'DUMPWriter', 'DUMPIOError',
           'DUMPFormatSpec']

#attr_map = {'id': 'atomID'}
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
    def __init__(self, fpath, **kwargs):
        super(DUMPReader, self).__init__(fpath=fpath, **kwargs)

        if self.fpath is not None:
            self.read()

    def read(self):
        """Read dump file."""
        pass
        #with open(self.fpath, 'r') as f:
        #    pass


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


class Snap(object):
    pass


class aselect(object):
    def __init__(self, data):
        self.data = data

    def all(self, *args):
        if len(args) == 0:
            for s in self.data.snaps:
                if not s.tselect:
                    continue
                for i in xrange(s.natoms):
                    s.aselect[i] = 1
                s.nselect = s.natoms
        else:
            s = self.data.snaps[self.data.findtime(args[0])]
            for i in xrange(s.natoms):
                s.aselect[i] = 1
            s.nselect = s.natoms


class tselect(object):
    def __init__(self, data):
        self.data = data

    def all(self):
        data = self.data
        for snap in data.snaps:
            snap.tselect = 1
        data.nselect = len(data.snaps)
        data.aselect.all()
        print('{}/{} snapshots selected'.format(data.nselect, data.nsnaps))

    def one(self, n):
        data = self.data
        for snap in data.snaps:
            snap.tselect = 0
        i = data.findtime(n)
        data.snaps[i].tselect = 1
        data.nselect = 1
        data.aselect.all()
        print('{}/{} snapshots selected'.format(data.nselect, data.nsnaps))

    def none(self):
        data = self.data
        for snap in data.snaps:
            snap.tselect = 0
        data.nselect = 0
        print('{}/{} snapshots selected'.format(data.nselect, data.nsnaps))

    def skip(self, n):
        data = self.data
        count = n - 1
        for snap in data.snaps:
            if not snap.tselect:
                continue
            count += 1
            if count == n:
                count = 0
                continue
            snap.tselect = 0
            data.nselect -= 1
        data.aselect.all()
        print('{}/{} snapshots selected'.format(data.nselect, data.nsnaps))


class DUMPData(object):
    """Class for reading and writing structure data in LAMMPS dump format.

    Parameters
    ----------
    *args : one or more dump files

    """
    def __init__(self, *args):
        #super(DUMPData, self).__init__()
        #self.atoms = Atoms()
        #self.dumpattrs = []
        self.snaps = []
        self.nsnaps = 0
        self.nselect = 0
        self.names = {}
        self.aselect = aselect(self)
        self.tselect = tselect(self)

        self.dumpfiles = []
        for fpath in args:
            self.dumpfiles.append(fpath)

        if len(self.dumpfiles) == 0:
            raise ValueError('No dump file specified.')

        #if len(args) == 1:
        #    self.increment = 0
        #    self.read_all()
        #else:
        #    self.increment = 1
        #    self.nextfile = 0
        #    self.eof = 0

        self.read_all()

    def read_all(self):
        """Read all snapshots from each dump file."""
        for dumpfile in self.dumpfiles:
            with open(dumpfile) as f:
                snap = self.read_snapshot(f)
                while snap is not None:
                    self.snaps.append(snap)
                    print(snap.time, end=' ')
                    sys.stdout.flush()
                    snap = self.read_snapshot(f)
        print()

        self.snaps.sort(key=attrgetter('time'))
        self.cull()
        self.nsnaps = len(self.snaps)

        print("read {:d} snapshots".format(self.nsnaps))

        self.tselect.all()
        if self.names:
            print('Dumped Atom attributes: {}'.format(self.names2str()))
        else:
            print('No dump column assignments')

        if 'x' not in self.names or 'y' not in self.names or \
                'z' not in self.names:
            print('dump scaling status unknown')
        elif self.nsnaps > 0:
            if self.scale_original == 1:
                self.unscale()
            elif self.scale_original == 0:
                print('dump is already unscaled')
            else:
                print('dump scaling status unknown')

    def read_snapshot(self, f):
        try:
            snap = Snap()
            f.readline()
            snap.time = int(f.readline().strip().split()[0])
            f.readline()
            snap.natoms = int(f.readline().strip())
            snap.aselect = np.zeros(snap.natoms)

            item = f.readline().strip()
            try:
                snap.boxstr = item.split('BOUNDS')[1].strip()
            except IndexError:
                snap.boxstr = ''

            snap.triclinic = False
            if 'xy' in snap.boxstr:
                snap.triclinic = True

            for axis, sf in zip(('x', 'y', 'z'), ('xy', 'xz', 'yz')):
                bounds = f.readline().strip().split()
                setattr(snap, axis + 'lo', float(bounds[0]))
                setattr(snap, axis + 'hi', float(bounds[1]))
                try:
                    setattr(snap, sf, float(bounds[-1]))
                except IndexError:
                    setattr(snap, sf, 0.0)

            if self.names:
                f.readline()
            else:
                self.scale_original = -1
                xflag = yflag = zflag = -1
                attrs = f.readline().strip().split()[2:]
                for i, attr in enumerate(attrs):
                    if attr in ('x', 'xu', 'xs', 'xsu'):
                        self.names['x'] = i
                        if attr in ('x', 'xu'):
                            xflag = 0
                        else:
                            xflag = 1
                    elif attr in ('y', 'yu', 'ys', 'ysu'):
                        self.names['y'] = i
                        if attr in ('y', 'yu'):
                            yflag = 0
                        else:
                            yflag = 1
                    elif attr in ('z', 'zu', 'zs', 'zsu'):
                        self.names['z'] = i
                        if attr in ('z', 'zu'):
                            zflag = 0
                        else:
                            zflag = 1
                    else:
                        self.names[attr] = i
                if xflag == yflag == zflag == 0:
                    self.scale_original = 0
                if xflag == yflag == zflag == 1:
                    self.scale_original = 1

                self.atomattrs = sorted(self.names, key=self.names.__getitem__)
                if 'id' in self.atomattrs:
                    self.atomattrs[self.atomattrs.index('id')] = 'atomID'

                if 'type' in self.atomattrs:
                    self.atomattrs[self.atomattrs.index('type')] = 'atomtype'

                self.attr_dtypes = [attr_dtypes[attr] if attr in attr_dtypes
                                    else float for attr in self.atomattrs]

            #atoms = Atoms()
            atoms = np.zeros((snap.natoms, len(self.atomattrs)), dtype=float)
            for n in xrange(snap.natoms):
                #line = [dtype(value) for dtype, value in
                #        zip(self.attr_dtypes, f.readline().strip().split())]
                #atoms.append(Atom(**dict(zip(self.atomattrs, line))))
                line = [float(attr) for attr in f.readline().strip().split()]
                atoms[n] = line

            snap.atoms = atoms
            return snap

        except Exception:
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
            x = self.names['x']
            y = self.names['y']
            z = self.names['z']
            for snap in self.snaps:
                self.scale_one(snap, x, y, z)

    def scale_one(self, snap, x, y, z):
        pass

    def unscale(self, ts=None):
        if ts is None:
            x = self.names['x']
            y = self.names['y']
            z = self.names['z']
            for snap in self.snaps:
                self.unscale_one(snap, x, y, z)

    def unscale_one(self, snap, x, y, z):
        if snap.xy == snap.xz == snap.yz == 0.0:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + atoms[:, x] * xprd
                atoms[:, y] = snap.ylo + atoms[:, y] * yprd
                atoms[:, z] = snap.zlo + atoms[:, z] * zprd
        else:
            xlo_bound = snap.xlo
            xhi_bound = snap.xhi
            ylo_bound = snap.ylo
            yhi_bound = snap.yhi
            zlo_bound = snap.zlo
            zhi_bound = snap.zhi
            xy = snap.xy
            xz = snap.xz
            yz = snap.yz
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
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + \
                    atoms[:, x] * h0 + atoms[:, y] * h5 + atoms[:, z] * h4
                atoms[:, y] = snap.ylo + atoms[:, y] * h1 + atoms[:, z] * h3
                atoms[:, z] = snap.zlo + atoms[:, z] * h2

    def wrap(self):
        pass

    def unwrap(self):
        pass

    def owrap(self):
        pass

    def atom(self, n, fields=None):
        pass

    def names2str(self):
        #return self.dumpattrs
        return ' '.join(sorted(self.names, key=self.names.__getitem__))

    def sort(self, *list):
        pass

    def sort_one(self, snap, id):
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

    def time(self):
        pass

    def vecs(self, n, *list):
        pass

    def newcolumn(self, str):
        pass

    def cull(self):
        i = 1
        while i < len(self.snaps):
            if self.snaps[i].time == self.snaps[i-1].time:
                del self.snaps[i]
            else:
                i += 1

    def delete(self):
        pass

    def extra(self):
        pass

    def findtime(self, n):
        pass

    def map(self, *pairs):
        pass

    def maxbox(self):
        pass

    def maxtype(self):
        pass

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
