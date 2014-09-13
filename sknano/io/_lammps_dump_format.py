# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS dump format (:mod:`sknano.io._lammps_dump_format`)
====================================================================

.. currentmodule:: sknano.io._lammps_dump_format

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import re

import numpy as np

from sknano.core import get_fpath

from ._base import Atom, StructureIO, StructureIOError, StructureFormatSpec, \
    default_comment_line

__all__ = ['DUMPData', 'DUMPReader', 'DUMPWriter', 'DUMPIOError',
           'DUMPFormatSpec']


class DUMPReader(StructureIO):
    """Class for reading `LAMMPS dump` file format.

    Parameters
    ----------
    fpath : str
        LAMMPS dump file path

    """
    def __init__(self, fpath=None, **kwargs):
        super(DUMPReader, self).__init__(fpath=fpath, **kwargs)

        if self.fpath is not None:
            self.read()

    def read(self):
        """Read dump file."""
        #with open(self.fpath, 'r') as f:
        #    pass

        self._parse_atoms()
        self._parse_atomtypes()
        self._parse_boxbounds()


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
                for i in xrange(s.Natoms):
                    s.aselect[i] = 1
                s.nselect = s.Natoms
        else:
            s = self.data.snaps[self.data.findtime(args[0])]
            for i in xrange(s.Natoms):
                s.aselect[i] = 1
            s.nselect = s.Natoms


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


class DUMPData(DUMPReader):
    """Class for reading and writing structure data in LAMMPS dump format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, files=None):
        super(DUMPData, self).__init__()
        self.snaps = []
        self.nsnaps = self.nselect = 0
        self.names = {}
        self.aselect = aselect(self)
        self.tselect = tselect(self)

    def atom(self, n, fields=None):
        pass

    def read_all(self):
        """Read all snapshots from each dump file."""
        pass

    def next(self):
        """Read next snapshot from list of dump files."""
        pass

    def read_snapshot(self, f):
        pass

    def scale(self, *list):
        pass

    def scale_one(self, snap, x, y, z):
        pass

    def unscale(self, *list):
        pass

    def unscale_one(self, snap, x, y, z):
        pass

    def wrap(self):
        pass

    def unwrap(self):
        pass

    def owrap(self):
        pass

    def names2str(self):
        pass

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

    def compare_atom(self, a, b):
        pass

    def compare_time(self, a, b):
        pass

    def cull(self):
        pass

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
