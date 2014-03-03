# -*- coding: utf-8 -*-
"""
======================================================
JSON format (:mod:`sknano.structure_io._json_format`)
======================================================

.. currentmodule:: sknano.structure_io._json_format

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

from ..chemistry import Atom, Atoms
from ..tools import get_fpath

from ._structure_data import StructureReader, StructureReaderError, \
    StructureWriter, default_comment_line

__all__ = ['JSONDATA', 'JSONReader', 'JSONWriter']


class JSONReader(StructureReader):
    """Class for reading json file format.

    Parameters
    ----------
    jsonfile : str
        json file

    """
    def __init__(self, fname=None):
        super(JSONReader, self).__init__(fname=fname)

        if fname is not None:
            self.read()

    def read(self):
        with open(self._fname, 'r') as f:
            Natoms = int(f.readline().strip())
            self._comment_line = f.readline().strip()
            lines = f.readlines()
            for line in lines:
                s = line.strip().split()
                if len(s) != 0:
                    atom = \
                        Atom(s[0], x=float(s[1]), y=float(s[2]), z=float(s[3]))
                    self._structure_atoms.append(atom)
            if self._structure_atoms.Natoms != Natoms:
                error_msg = '`jsonfile` contained {} atoms '.format(
                    self._structure_atoms.Natoms) + 'but should contain ' + \
                    '{}'.format(Natoms)
                raise StructureReaderError(error_msg)


class JSONWriter(StructureWriter):
    """Class for writing json chemical file format."""

    @classmethod
    def write(cls, fname=None, atoms=None, comment_line=None):
        """Write structure data to file.

        Parameters
        ----------
        fname : str
        atoms : :py:class:`~sknano.chemistry.Atoms`
            An :py:class:`~sknano.chemistry.Atoms` instance.
        comment_line : str, optional
            A string written to the first line of `json` file. If `None`,
            then it is set to the full path of the output `json` file.

        """
        if not isinstance(atoms, Atoms):
            raise TypeError('atoms argument must be an `Atoms` instance')
        else:
            fname = get_fpath(fname=fname, ext='json', overwrite=True,
                              add_fnum=False)
            if comment_line is None:
                comment_line = default_comment_line

            atoms.rezero_coords()

            with open(fname, 'w') as f:
                f.write('{:d}\n'.format(atoms.Natoms))
                f.write('{}\n'.format(comment_line))
                for atom in atoms:
                    f.write('{:>3s}{:15.8f}{:15.8f}{:15.8f}\n'.format(
                        atom.symbol, atom.x, atom.y, atom.z))


class JSONDATA(JSONReader):
    """Class for reading and writing structure data in JSON data format.

    Parameters
    ----------
    fname : str, optional

    """
    def __init__(self, fname=None):
        super(JSONDATA, self).__init__(fname=fname)

    def write(self, jsonfile=None):
        """Write json file.

        Parameters
        ----------
        jsonfile : {None, str}, optional

        """
        try:
            if (jsonfile is None or jsonfile == '') and \
                    (self.fname is None or self.fname == ''):
                error_msg = '`jsonfile` must be a string at least 1 ' + \
                    'character long.'
                if jsonfile is None:
                    raise TypeError(error_msg)
                else:
                    raise ValueError(error_msg)
            else:
                jsonfile = self._fname
            JSONWriter.write(fname=jsonfile, atoms=self._structure_atoms,
                             comment_line=self._comment_line)
        except (TypeError, ValueError) as e:
            print(e)
