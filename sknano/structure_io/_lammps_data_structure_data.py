# -*- coding: utf-8 -*-
"""
==============================================================================
LAMMPS data format (:mod:`sknano.structure_io._lammps_data_structure_data`)
==============================================================================

.. currentmodule:: sknano.structure_io._lammps_data_structure_data

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

import numpy as np

from pksci.chemistry import Atom
from pkshared.tools.fiofuncs import get_fpath

from ._structure_data import StructureReader, StructureWriter


__all__ = ['DATAReader', 'DATAWriter', 'LAMMPSDATA']


class DATAReader(StructureReader):
    """Class for reading ``LAMMPS data`` file format.

    Parameters
    ----------
    datafile : str
        LAMMPS data file
    atom_style : {'full', 'atomic'}, optional

    """
    def __init__(self, datafile=None, atom_style='full'):
        super(DATAReader, self).__init__(fname=datafile)

        from ._structure_specs import LAMMPSDATASpecs
        data_specs = LAMMPSDATASpecs(atom_style=atom_style)
        self._data_headers = data_specs.properties['headers']
        self._data_sections = data_specs.properties['sections']
        self._section_properties = data_specs.section_properties
        self._section_syntax_dict = data_specs.section_syntax_dict

        self._headers = {}
        self._sections = {}
        self._boxbounds = {}

        if datafile is not None:
            self._read()
            self._parse_atoms()
            self._Natoms = self._atoms.Natoms
            self._parse_boxbounds()

    @property
    def headers(self):
        """DATA file headers."""
        return self._headers

    @property
    def sections(self):
        """DATA file sections."""
        return self._sections

    @property
    def boxbounds(self):
        """Box bounds."""
        return self._boxbounds

    def _read(self):
        """Read data file."""
        with open(self._fname, 'r') as f:
            self._comment_line = f.readline().strip()

            while True:
                line = f.readline().strip()
                if len(line) == 0:
                    continue
                found = False
                for key in self._data_headers.iterkeys():
                    if key in line:
                        found = True
                        self._headers[key] = \
                            [self._data_headers[key]['dtype'](s) for s in
                                [[ss for ss in line.split()][i] for i in
                                 range(self._data_headers[key]['items'])]]
                        if len(self._headers[key]) == 1:
                            self._headers[key] = self._headers[key][0]
                        break
                if not found:
                    break

            while True:
                found = False
                for section_key, header_key in self._data_sections.iteritems():
                    if section_key in line:
                        found = True
                        f.readline()
                        Nitems = self._headers[header_key]
                        data = []
                        for n in xrange(Nitems):
                            tmp = []
                            line = f.readline().strip().split()
                            for i, props in \
                                enumerate(self._section_properties[
                                    section_key].itervalues()):
                                tmp.append(props['dtype'](line[i]))
                            data.append(tmp)
                        self._sections[section_key] = data[:]
                        #self._sections[section_key] = \
                        #    [[props['dtype'](s) for props in
                        #        self._section_properties[
                        #            section_key].itervalues() for
                        #        s in f.readline().strip().split()] for n in
                        #        xrange(Nitems)]
                        break
                f.readline()
                line = f.readline().strip()
                if len(line) == 0:
                    break

    def _parse_atoms(self):
        """Populate Atoms object with Atom objects"""
        atoms = self._sections['Atoms']
        masses = self._sections['Masses']
        velocities = self._sections['Velocities']
        atom_kwargs = {'atomID': None, 'moleculeID': None,
                       'q': None, 'atomtype': None, 'mass': None,
                       'x': None, 'y': None, 'z': None,
                       'vx': None, 'vy': None, 'vz': None}
        atoms_section_syntax = self._section_syntax_dict['Atoms']
        masses_section_syntax = self._section_syntax_dict['Masses']
        velocities_section_syntax = self._section_syntax_dict['Velocities']

        for atom in atoms:
            for kw in atom_kwargs.iterkeys():
                if kw in atoms_section_syntax:
                    atom_kwargs[kw] = \
                        atom[self._section_properties['Atoms'][kw]['index']]
                elif kw in masses_section_syntax:
                    atomtype = \
                        atom[self._section_properties[
                            'Atoms']['atomtype']['index']]
                    atom_kwargs[kw] = \
                        masses[atomtype-1][self._section_properties[
                            'Masses'][kw]['index']]
                elif kw in velocities_section_syntax and \
                        len(velocities) == len(atoms):
                    atomID = \
                        atom[self._section_properties[
                            'Atoms']['atomID']['index']]
                    for velocity in velocities:
                        velocity_atomID = \
                            velocity[self._section_properties[
                                'Velocities']['atomID']['index']]
                        if velocity_atomID == atomID:
                            atom_kwargs[kw] = \
                                velocity[self._section_properties[
                                    'Velocities'][kw]['index']]
                else:
                    print('unknown atom keyword: {}'.format(kw))

            _atom = Atom(**atom_kwargs)
            self._atoms.append(_atom)

    #def _parse_atomtypes(self):
    #    mass_syntax = self._section_properties['Masses']['atomtype']

    def _parse_boxbounds(self):
        for dim in ('x', 'y', 'z'):
            bounds = \
                self._headers[' '.join([dim + lim for lim in ('lo', 'hi')])]
            self._boxbounds[dim] = {'min': bounds[0], 'max': bounds[-1]}

    def get(self, section_key, colnum=None, colname=None, colindex=None):
        """Return section with ``section key``.

        Parameters
        ----------
        section_key : str
        colnum : int, optional
        colname : str, optional
        colindex : int, optional

        Returns
        -------
        sequence

        """
        section_data = None
        try:
            section_data = self._sections[section_key]
            section_syntax = self._section_syntax_dict[section_key]
        except KeyError as e:
            print(e)
        else:
            try:
                colidx = None
                if colnum is not None:
                    colidx = int(colnum - 1)
                elif colname is not None:
                    colidx = \
                        self._section_properties[section_key][colname]['index']
                elif colindex is not None:
                    colidx = int(colindex)
            except (KeyError, TypeError, ValueError) as e:
                print(e)
            else:
                try:
                    colname = section_syntax[colidx]
                    coltype = \
                        self._section_properties[section_key][colname]['dtype']
                    section_data = \
                        np.asarray(
                            section_data, dtype=coltype)[:, colidx].tolist()
                except Exception as e:
                    print(e)
        finally:
            return section_data


class DATAWriter(StructureWriter):
    """Class for writing LAMMPS data chemical file format."""

    @classmethod
    def write(cls, fname=None, atoms=None, boxbounds=None, comment_line=None):
        """Write structure data to file.

        Parameters
        ----------
        fname : str
        atoms : Atoms
            :py:class:`~pksci.chemistry.Atoms` instance.
        boxbounds : dict, optional
            If ``None``, determined automatically from atom coordinates.
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


class LAMMPSDATA(DATAReader):
    """Class for reading and writing structure data in LAMMPS data format.

    Parameters
    ----------
    datafile : str, optional

    """
    def __init__(self, datafile=None):
        super(LAMMPSDATA, self).__init__(datafile=datafile)

    def map_colinfo(self):
        pass

    def maxbox(self):
        pass

    def maxtype(self):
        pass

    def replace(self, section_key):
        pass

    def write(self):
        DATAWriter.write(fname=self._fname, atoms=self._atoms,
                         boxbounds=self._boxbounds,
                         comment_line=self._comment_line)

    def newxyz(self):
        pass

    def delete(self):
        pass
