# -*- coding: utf-8 -*-
"""
===========================================================================
Structure readers (:mod:`sknano.structure_io._writers`)
===========================================================================

.. currentmodule:: sknano.structure_io._writers

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from abc import ABCMeta, abstractmethod

from pkshared.tools.fiofuncs import get_fpath

__all__ = ['StructureWriter', 'DATAWriter', 'XYZWriter']


class StructureWriter(object):
    __metaclass__ = ABCMeta
    """Abstract superclass for writing structure data.

    Parameters
    ----------
    fname : str
        structure file

    """
    def __init__(self, fname=None):
        self._fname = fname

    @property
    def fname(self):
        return self._fname

    @abstractmethod
    def write(self):
        """Read in structure data from file"""
        return NotImplemented


class DATAWriter(StructureWriter):
    """Class for writing LAMMPS data chemical file format.

    Parameters
    ----------
    datafile : str
        LAMMPS data file

    """
    def __init__(self, datafile):
        super(DATAWriter, self).__init__(fname=datafile)

    @classmethod
    def write(cls, fname=None, atoms=None, boxbounds=None, comment_line=None):
        """Write structure data to file.

        Parameters
        ----------
        fname : str
        atoms : Atoms
        boxbounds : dict, optional
            ..versionchanged:: 0.3.32
            No longer required. If ``None``, determined automatically from
            atom coordinates.
        comment_line : str, optional

        """
        if fname is None:
            raise TypeError('fname argument must be a string!')
        elif atoms is None:
            raise TypeError('atoms argument must be an Atoms object')
        else:
            fname = get_fpath(fname=fname, ext='data', overwrite=True)
            if comment_line is None:
                comment_line = fname

            atoms.fix_minus_zero_coords()
            atomtypes = atoms.atomtypes

            Natoms = atoms.Natoms
            Natoms_width = \
                8 if len(str(Natoms)) <= 12 else len(str(Natoms)) + 4
            Ntypes = atoms.Ntypes
            Ntypes_width = Natoms_width

            if boxbounds is None:
                boxbounds = {'x': {'min': None, 'max': None},
                             'y': {'min': None, 'max': None},
                             'z': {'min': None, 'max': None}}

                for i, dim in enumerate(('x', 'y', 'z')):
                    boxbounds[dim]['min'] = atoms.coords[:, i].min()
                    boxbounds[dim]['max'] = atoms.coords[:, i].max()

            lohi_width = 0
            for dim in ('x', 'y', 'z'):
                lohi_width = \
                    max(lohi_width, len('{:.6f} {:.6f}'.format(
                        boxbounds[dim]['min'], boxbounds[dim]['max'])) + 4)

            with open(fname, 'w') as f:
                f.write('# {}\n\n'.format(comment_line))
                f.write('{}atoms\n'.format(
                    '{:d}'.format(Natoms).ljust(Natoms_width)))
                f.write('{}atom types\n'.format(
                    '{:d}'.format(Ntypes).ljust(Ntypes_width)))
                f.write('\n')
                for dim in ('x', 'y', 'z'):
                    f.write('{}{dim}lo {dim}hi\n'.format(
                        '{:.6f} {:.6f}'.format(
                            boxbounds[dim]['min'],
                            boxbounds[dim]['max']).ljust(lohi_width),
                        dim=dim))
                f.write('\nMasses\n\n')
                for atomtype, properties in atomtypes.iteritems():
                    f.write('{}{:.4f}\n'.format(
                        '{:d}'.format(atomtype).ljust(Natoms_width),
                        properties['mass']))
                f.write('\nAtoms\n\n')
                for atomID, atom in enumerate(atoms, start=1):
                    line = ''
                    line += "{:>{}}".format(atomID, len(str(Natoms)) + 1)
                    line += "{:>{}}".format(atom.moleculeID, 3)
                    line += "{:>{}}".format(atom.atomtype,
                                            len(str(Ntypes)) + 1)
                    line += "{:>{}}".format('{:.1f}'.format(atom.q), 4)
                    line += "{:>{}}".format('{:f}'.format(atom.x), 14)
                    line += "{:>{}}".format('{:f}'.format(atom.y), 14)
                    line += "{:>{}}".format('{:f}'.format(atom.z), 14)
                    line += "{:>{}}".format('{:d}'.format(atom.nx), 3)
                    line += "{:>{}}".format('{:d}'.format(atom.ny), 3)
                    line += "{:>{}}\n".format('{:d}'.format(atom.nz), 3)

                    f.write(line)

                f.write('\nVelocities\n\n')
                for atomID, atom in enumerate(atoms, start=1):
                    line = ''
                    line += "{:>{}}".format(atomID, len(str(Natoms)) + 1)
                    line += "{:>{}}".format('{:f}'.format(atom.vx), 14)
                    line += "{:>{}}".format('{:f}'.format(atom.vy), 14)
                    line += "{:>{}}\n".format('{:f}'.format(atom.vz), 14)

                    f.write(line)


class XYZWriter(StructureWriter):
    """Class for writing xyz chemical file format.

    Parameters
    ----------
    xyzfile : str
        xyz structure file

    """
    def __init__(self, xyzfile):
        super(XYZWriter, self).__init__(fname=xyzfile)

    @classmethod
    def write(cls, fname=None, atoms=None, comment_line=None):
        """Write structure data to file.

        Parameters
        ----------
        fname : str
        atoms : Atoms
        comment_line : str, optional

        """
        if fname is None:
            raise TypeError('fname argument must be a string!')
        elif atoms is None:
            raise TypeError('atoms argument must be an Atoms object')
        else:
            fname = get_fpath(fname=fname, ext='xyz', overwrite=True)
            if comment_line is None:
                comment_line = fname

            with open(fname, 'w') as f:
                f.write('{:d}\n'.format(atoms.Natoms))
                f.write('{}\n'.format(comment_line))
                for atom in atoms.atomlist:
                    atom.fix_minus_zero_coords()
                    f.write('{:3s} {:10.5f} {:10.5f} {:10.5f}\n'.format(
                        atom.symbol, atom.x, atom.y, atom.z))
