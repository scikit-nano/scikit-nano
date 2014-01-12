# -*- coding: utf-8 -*-
"""
==============================================================================
XYZ format (:mod:`sknano.structure_io._xyz_structure_data`)
==============================================================================

.. currentmodule:: sknano.structure_io._xyz_structure_data

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from pkshared.tools.fiofuncs import get_fpath

from ..chemistry import Atom, Atoms
from ._structure_data import StructureReader, StructureReaderError, \
    StructureWriter

__all__ = ['XYZDATA', 'XYZReader', 'XYZWriter']


class XYZReader(StructureReader):
    """Class for reading xyz chemical file format.

    Parameters
    ----------
    xyzfile : str
        xyz structure file

    """
    def __init__(self, fname=None):
        super(XYZReader, self).__init__(fname=fname)

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
                    self._atoms.append(atom)
            if self._atoms.Natoms != Natoms:
                error_msg = '`xyzfile` contained {} atoms '.format(
                    self._atoms.Natoms) + 'but should contain ' + \
                    '{}'.format(Natoms)
                raise StructureReaderError(error_msg)


class XYZWriter(StructureWriter):
    """Class for writing xyz chemical file format."""

    @classmethod
    def write(cls, fname=None, atoms=None, comment_line=None):
        """Write structure data to file.

        Parameters
        ----------
        fname : str
        atoms : `Atoms`
            An :py:class:`Atoms` instance.
        comment_line : str, optional
            A string written to the first line of ``xyz`` file. If ``None``,
            then it is set to the full path of the output ``xyz`` file.

        """
        if not isinstance(atoms, Atoms):
            raise TypeError('atoms argument must be an ``Atoms`` instance')
        else:
            fname = get_fpath(fname=fname, ext='xyz', overwrite=True,
                              add_fnum=False)
            if comment_line is None:
                comment_line = fname

            atoms.fix_minus_zero_coords()

            with open(fname, 'w') as f:
                f.write('{:d}\n'.format(atoms.Natoms))
                f.write('{}\n'.format(comment_line))
                for atom in atoms:
                    f.write('{:>3s}{:15.8f}{:15.8f}{:15.8f}\n'.format(
                        atom.symbol, atom.x, atom.y, atom.z))


class XYZDATA(XYZReader):
    """Class for reading and writing structure data in XYZ data format.

    Parameters
    ----------
    fname : str, optional

    """
    def __init__(self, fname=None):
        super(XYZDATA, self).__init__(fname=fname)

    def write(self, xyzfile=None):
        """Write xyz file.

        Parameters
        ----------
        xyzfile : {None, str}, optional

        """
        try:
            if (xyzfile is None or xyzfile == '') and \
                    (self.fname is None or self.fname == ''):
                error_msg = '`xyzfile` must be a string at least 1 ' + \
                    'character long.'
                if xyzfile is None:
                    raise TypeError(error_msg)
                else:
                    raise ValueError(error_msg)
            else:
                xyzfile = self._fname
            XYZWriter.write(fname=xyzfile, atoms=self._atoms,
                            comment_line=self._comment_line)
        except (TypeError, ValueError) as e:
            print(e)
