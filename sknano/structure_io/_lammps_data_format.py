# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS data format (:mod:`sknano.structure_io._lammps_data_format`)
====================================================================

.. currentmodule:: sknano.structure_io._lammps_data_format

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
import os

import numpy as np

from .atoms import LAMMPSAtom as Atom, LAMMPSAtoms, lammps_atom_styles
from ..tools import get_fpath

from ._structure_data import StructureReader, StructureWriter, \
    StructureConverter, StructureFormat, StructureDataError, \
    default_comment_line

__all__ = ['DATAData', 'DATAReader', 'DATAWriter', 'DATA2XYZConverter',
           'DATAFormat', 'DATAError', 'LAMMPSDATA']


class DATAReader(StructureReader):
    """`StructureReader` class for reading `LAMMPS data` file format.

    Parameters
    ----------
    fpath : str
        `LAMMPS data` file path
    atom_style : {'full', 'atomic'}, optional

    """
    def __init__(self, fpath=None, atom_style='full'):
        super(DATAReader, self).__init__(fpath=fpath)

        self._structure_atoms = LAMMPSAtoms()

        data_format = DATAFormat(atom_style=atom_style)
        self._data_headers = data_format.properties['headers']
        self._data_sections = data_format.properties['sections']
        self._section_properties = data_format.section_properties
        self._section_syntax_dict = data_format.section_syntax_dict

        self._headers = {}
        self._sections = {}
        self._boxbounds = {}

        if fpath is not None:
            self.read()

    @property
    def headers(self):
        """LAMMPS data file headers."""
        return self._headers

    @property
    def sections(self):
        """LAMMPS data file sections."""
        return self._sections

    @property
    def boxbounds(self):
        """Box bounds."""
        return self._boxbounds

    def read(self):
        """Read data file."""
        with open(self.fpath, 'r') as f:
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
                            [self._data_headers[key]['dtype'](float(s)) for s
                             in [[ss for ss in line.split()][i] for i in
                                 range(self._data_headers[key]['items'])]]
                        if len(self._headers[key]) == 1:
                            # since the list contains only one element,
                            # replace list value with element 0
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
                        #self._sections[section_key] = \
                        #    [[props['dtype'](s) for props in
                        #        self._section_properties[
                        #            section_key].itervalues() for
                        #        s in f.readline().strip().split()] for n in
                        #        xrange(Nitems)]
                        data = []
                        for n in xrange(Nitems):
                            tmp = []
                            line = f.readline().strip().split()
                            for i, props in enumerate(
                                self._section_properties[
                                    section_key].itervalues()):
                                try:
                                    tmp.append(props['dtype'](float(line[i])))
                                except IndexError:
                                    break
                            data.append(tmp)
                        self._sections[section_key] = data[:]
                        break
                f.readline()
                line = f.readline().strip()
                if len(line) == 0:
                    break

        self._parse_atoms()
        self._parse_atomtypes()
        self._parse_boxbounds()

    def _parse_atoms(self):
        """Populate `Atoms` object with `Atom` objects"""
        atoms_section = self._sections['Atoms']
        atoms_section_syntax = self._section_syntax_dict['Atoms']

        masses_section = self._sections['Masses']
        masses_section_syntax = self._section_syntax_dict['Masses']

        try:
            velocities_section = self._sections['Velocities']
            velocities_section_syntax = self._section_syntax_dict['Velocities']
        except KeyError:
            velocities_section = []
            velocities_section_syntax = {}

        atom_kwargs = {'atomID': None, 'moleculeID': None,
                       'q': None, 'atomtype': None, 'mass': None,
                       'x': None, 'y': None, 'z': None,
                       'vx': None, 'vy': None, 'vz': None}

        for lmps_atom in atoms_section:
            for kw in atom_kwargs.iterkeys():
                if kw in atoms_section_syntax:
                    atom_kwargs[kw] = \
                        lmps_atom[
                            self._section_properties['Atoms'][kw]['index']]
                elif kw in masses_section_syntax:
                    atomtype = \
                        lmps_atom[
                            self._section_properties[
                                'Atoms']['atomtype']['index']]
                    atom_kwargs[kw] = \
                        masses_section[atomtype-1][
                            self._section_properties['Masses'][kw]['index']]
                elif kw in velocities_section_syntax and \
                        len(velocities_section) == len(atoms_section):
                    atomID = \
                        lmps_atom[
                            self._section_properties[
                                'Atoms']['atomID']['index']]
                    for velocity in velocities_section:
                        velocity_atomID = \
                            velocity[self._section_properties[
                                'Velocities']['atomID']['index']]
                        if velocity_atomID == atomID:
                            atom_kwargs[kw] = \
                                velocity[self._section_properties[
                                    'Velocities'][kw]['index']]
                #else:
                #    print('unknown atom keyword: {}'.format(kw))

            atom = Atom(**atom_kwargs)
            self._structure_atoms.append(atom)

    def _parse_atomtypes(self):
        Ntypes = self._structure_atoms.Ntypes
        atomtypes = self._structure_atoms.atomtypes
        if Ntypes != self._headers['atom types']:
            for atomtype in xrange(1, self._headers['atom types'] + 1):
                if atomtype not in atomtypes:
                    mass = self._sections['Masses'][atomtype - 1][
                        self._section_properties['Masses']['mass']['index']]
                    self._structure_atoms.add_atomtype(
                        Atom(atomtype=atomtype, mass=mass))

    def _parse_boxbounds(self):
        for dim in ('x', 'y', 'z'):
            bounds = \
                self._headers[' '.join([dim + lim for lim in ('lo', 'hi')])]
            self._boxbounds[dim] = {'min': bounds[0], 'max': bounds[-1]}

    def get(self, section_key, colnum=None, colname=None, colindex=None):
        """Return section with `section_key`.

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
                except TypeError:
                    pass
        finally:
            return section_data


class DATAWriter(StructureWriter):
    """`StructureWriter` class for writing `LAMMPS data` file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, atoms=None,
              atom_style='full', boxbounds=None, comment_line=None,
              assume_unique_atoms=False, enforce_consecutive_atomIDs=True,
              pad_box=True, xpad=10., ypad=10., zpad=10., verbose=False):
        """Write structure data to file.

        Parameters
        ----------
        fname : str, optional
            Output file name.
        outpath : str, optional
            Output file path.
        fpath : str, optional
            Full path (directory path + file name) to output data file.
        atoms : :class:`~sknano.chemistry.Atoms`
            An :class:`~sknano.chemistry.Atoms` instance.
        boxbounds : dict, optional
            If `None`, determined automatically from the `atoms` coordinates.
        comment_line : str, optional
            A string written to the first line of `data` file. If `None`,
            then it is set to the full path of the output `data` file.
        assume_unique_atoms : bool, optional
            Check that each :class:`~sknano.chemistry.Atom` in `atoms`
            has a unique :attr:`~sknano.chemistry.Atom.atomID`.
            If the check fails, then assign a unique
            :attr:`~sknano.chemistry.Atom.atomID` to each
            :class:`~sknano.chemistry.Atom`.
            If `assume_unique_atoms` is True, but the atomID's are not unique,
            LAMMPS will not be able to read the data file.
        verbose : bool, optional
            verbose output

        """
        if fpath is None:
            fpath = get_fpath(fname=fname, ext='data', outpath=outpath,
                              overwrite=True, add_fnum=False)
        if comment_line is None:
            comment_line = default_comment_line

        atoms.rezero_coords()

        atomtypes = atoms.atomtypes

        Natoms = atoms.Natoms
        Natoms_width = \
            8 if len(str(Natoms)) <= 12 else len(str(Natoms)) + 4
        Ntypes = atoms.Ntypes
        Ntypes_width = Natoms_width

        atomID_width = len(str(Natoms)) + 1
        atomtype_width = len(str(Ntypes)) + 1

        if (enforce_consecutive_atomIDs and
            atoms.atom_ids.max() != atoms.Natoms) or \
                (not assume_unique_atoms and
                 len(set(atoms.atom_ids)) != atoms.Natoms):
            atoms.assign_unique_ids()

        if boxbounds is None:
            boxbounds = {'x': {'min': None, 'max': None},
                         'y': {'min': None, 'max': None},
                         'z': {'min': None, 'max': None}}

            for i, dim in enumerate(('x', 'y', 'z')):
                boxbounds[dim]['min'] = atoms.coords[:, i].min()
                boxbounds[dim]['max'] = atoms.coords[:, i].max()

        boxpad = {'x': xpad, 'y': ypad, 'z': zpad}
        pad_eps = 0.001
        if pad_box:
            #for dim, pad in boxpad.iteritems():
            for i, dim in enumerate(('x', 'y', 'z')):
                pad = boxpad[dim]
                if abs(boxbounds[dim]['min'] - atoms.coords[:, i].min()) \
                        < pad - pad_eps:
                    boxbounds[dim]['min'] = boxbounds[dim]['min'] - pad
                if abs(boxbounds[dim]['max'] - atoms.coords[:, i].max()) \
                        < pad - pad_eps:
                    boxbounds[dim]['max'] = boxbounds[dim]['max'] + pad

        lohi_width = 0
        for dim in ('x', 'y', 'z'):
            lohi_width = \
                max(lohi_width, len('{:.6f} {:.6f}'.format(
                    boxbounds[dim]['min'], boxbounds[dim]['max'])) + 4)

        with open(fpath, 'w') as f:
            f.write('# {}\n\n'.format(comment_line.lstrip('#').strip()))

            f.write('{}atoms\n'.format(
                '{:d}'.format(Natoms).ljust(Natoms_width)))
            f.write('{}atom types\n\n'.format(
                '{:d}'.format(Ntypes).ljust(Ntypes_width)))

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
            for atom in atoms:
                line = ''
                line += "{:>{}}".format(atom.atomID, atomID_width)
                line += "{:>{}}".format(atom.moleculeID, 3)
                line += "{:>{}}".format(atom.atomtype, atomtype_width)
                line += "{:>{}}".format('{:.1f}'.format(atom.q), 4)
                line += "{:>{}}".format('{:f}'.format(atom.x), 14)
                line += "{:>{}}".format('{:f}'.format(atom.y), 14)
                line += "{:>{}}".format('{:f}'.format(atom.z), 14)
                line += "{:>{}}".format('{:d}'.format(atom.nx), 3)
                line += "{:>{}}".format('{:d}'.format(atom.ny), 3)
                line += "{:>{}}".format('{:d}'.format(atom.nz), 3)
                line += '\n'

                f.write(line)

            f.write('\nVelocities\n\n')
            for atom in atoms:
                line = ''
                line += "{:>{}}".format(atom.atomID, atomID_width)
                line += "{:>{}}".format('{:f}'.format(atom.vx), 14)
                line += "{:>{}}".format('{:f}'.format(atom.vy), 14)
                line += "{:>{}}".format('{:f}'.format(atom.vz), 14)
                line += '\n'

                f.write(line)


class DATA2XYZConverter(StructureConverter):
    """
    `StructureConverter` class for converting `LAMMPS data` to `xyz` format.

    .. versionadded:: 0.2.9

    Parameters
    ----------
    datafile : str

    """
    def __init__(self, datafile):
        self._datafile = datafile
        self._xyzfile = os.path.splitext(self._datafile)[0] + '.xyz'

        super(DATA2XYZConverter, self).__init__(
            infile=self._datafile, outfile=self._xyzfile)

    @property
    def datafile(self):
        """`LAMMPS data` file."""
        return self.infile

    @property
    def xyzfile(self):
        """`xyz` file name."""
        return self.outfile

    def convert(self, return_reader=False):
        """Convert `LAMMPS data` to `xyz` chemical file format.

        Parameters
        ----------
        return_reader : bool, optional
            Return an instance of `XYZReader`

        Returns
        -------
        `XYZReader` (only if `return_reader` is True)

        """
        from .atoms import AtomsConverter
        from ._xyz_format import XYZReader, XYZWriter

        datareader = DATAReader(fpath=self.infile)
        atoms = AtomsConverter(atoms=datareader.atoms, to='xyz').atoms
        comment_line = datareader.comment_line

        XYZWriter.write(fpath=self.outfile, atoms=atoms,
                        comment_line=comment_line)

        if return_reader:
            return XYZReader(fpath=self.outfile)


class DATAData(DATAReader):
    """Class for reading and writing `StructureData` in `LAMMPS data` format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None):
        super(DATAData, self).__init__(fpath=fpath)

    def delete(self, key):
        if key in self.headers:
            del self._headers[key]
        elif key in self.sections:
            del self._sections[key]
        else:
            raise DATAError('{} not in DATA object'.format(key))

    def map(self, *pairs):
        pass

    def maxbox(self):
        pass

    def maxtype(self):
        pass

    def newxyz(self):
        pass

    def reorder(self, colname, *order):
        pass

    def replace(self, section_key, new_data, colnum=None,
                colname=None, colindex=None):
        """Replace section data.

        Parameters
        ----------
        section_key : str
        new_data : sequence
        colnum : int, optional
        colname : str, optional
        colindex : int, optional

        """
        colidx = None

        # for backwards compatibility with the pizza.py data module,
        # first check positional arguments to see if this method was called
        # using the pizza.py data module signature which expects positional
        # arguments of the following type:
        #   data.replace(str, int, list)
        if isinstance(new_data, (int, float)) and \
                isinstance(colnum, (np.ndarray, list)):
            colidx = int(new_data) - 1
            new_data = np.asarray(colnum)
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
                raise DATAError('replace called with invalid arguments')
        atom_attr = self._section_syntax_dict[section_key][colidx]
        attr_dtype = \
            self._section_properties[section_key][atom_attr]['dtype']
        new_data = np.asarray(new_data, dtype=attr_dtype)

        for i, atom in enumerate(self._structure_atoms):
            self._sections[section_key][i][colidx] = \
                attr_dtype(float(new_data[i]))
            setattr(atom, atom_attr, attr_dtype(float(new_data[i])))

    def viz(self, isnap):
        pass

    def write(self, datafile=None):
        """Write data file.

        Parameters
        ----------
        datafile : {None, str}, optional

        """
        try:
            if (datafile is None or datafile == '') and \
                    (self.fpath is None or self.fpath == ''):
                error_msg = \
                    '`datafile` must be a string at least 1 character long.'
                if datafile is None:
                    raise TypeError(error_msg)
                else:
                    raise ValueError(error_msg)
            else:
                datafile = self.fpath

            DATAWriter.write(fname=datafile,
                             atoms=self._structure_atoms,
                             boxbounds=self._boxbounds,
                             comment_line=self._comment_line)
        except (TypeError, ValueError) as e:
            print(e)

    @classmethod
    def format_spec(cls, atom_style='full'):
        return DATAFormat(atom_style=atom_style)

LAMMPSDATA = DATAData


class DATAError(StructureDataError):
    """`StructureDataError` `Exception` for `LAMMPS data` file format.

    Parameters
    ----------
    msg : str
        Error message

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class DATAFormat(StructureFormat):
    """`StructureFormat` class defining properties for `LAMMPS data` format.

    Parameters
    ----------
    atom_style : {'full'}, optional
        LAMMPS atom style.

    """
    def __init__(self, atom_style='full'):
        super(DATAFormat, self).__init__()
        self._atom_style = atom_style
        self._section_properties = OrderedDict()

        atoms_section_syntax = {}
        for atom_style, var_list in lammps_atom_styles.iteritems():
            atoms_section_syntax[atom_style] = []
            for var in var_list:
                atoms_section_syntax[atom_style].append(var.replace('-', ''))

        velocities_section_syntax = {}
        velocities_section_syntax['full'] = ['atomID', 'vx', 'vy', 'vz']

        syntax_dtypes = \
            {'atomID': int, 'moleculeID': int, 'atomtype': int, 'q': float,
             'mass': float, 'x': float, 'y': float, 'z': float,
             'nx': int, 'ny': int, 'nz': int,
             'vx': float, 'vy': float, 'vz': float}

        section_keys = ['Atoms', 'Masses', 'Velocities']
        self._section_syntax_dict = \
            {'Atoms': atoms_section_syntax[self.atom_style],
             'Masses': ['atomtype', 'mass'],
             'Velocities': velocities_section_syntax[self.atom_style]}

        for section_key in section_keys:
            self._section_properties[section_key] = OrderedDict()
            section_syntax_list = self._section_syntax_dict[section_key]
            for i, syntax in enumerate(section_syntax_list):
                self._section_properties[section_key][syntax] = \
                    {'dtype': syntax_dtypes[syntax],
                     'colnum': i+1,
                     'index': i}

        self._header_keys = ['atoms',
                             'bonds',
                             'angles',
                             'dihedrals',
                             'impropers',
                             'atom types',
                             'bond types',
                             'angle types',
                             'dihedral types',
                             'improper types',
                             'extra bond per atom',
                             'extra angle per atom',
                             'extra dihedral per atom',
                             'extra improper per atom',
                             'ellipsoids',
                             'lines',
                             'triangles',
                             'bodies',
                             'xlo xhi',
                             'ylo yhi',
                             'zlo zhi',
                             'xy xz yz']

        self._headers = {'atoms': {'dtype': int, 'items': 1},
                         'bonds': {'dtype': int, 'items': 1},
                         'angles': {'dtype': int, 'items': 1},
                         'dihedrals': {'dtype': int, 'items': 1},
                         'impropers': {'dtype': int, 'items': 1},
                         'atom types': {'dtype': int, 'items': 1},
                         'bond types': {'dtype': int, 'items': 1},
                         'angle types': {'dtype': int, 'items': 1},
                         'dihedral types': {'dtype': int, 'items': 1},
                         'improper types': {'dtype': int, 'items': 1},
                         'extra bond per atom': {'dtype': int, 'items': 1},
                         'extra angle per atom': {'dtype': int, 'items': 1},
                         'extra dihedral per atom': {'dtype': int, 'items': 1},
                         'extra improper per atom': {'dtype': int, 'items': 1},
                         'ellipsoids': {'dtype': int, 'items': 1},
                         'lines': {'dtype': int, 'items': 1},
                         'triangles': {'dtype': int, 'items': 1},
                         'bodies': {'dtype': int, 'items': 1},
                         'xlo xhi': {'dtype': float, 'items': 2},
                         'ylo yhi': {'dtype': float, 'items': 2},
                         'zlo zhi': {'dtype': float, 'items': 2},
                         'xy xz yz': {'dtype': float, 'items': 3}}

        self._properties['headers'] = self._headers

        # A LAMMPS data file is partitioned into sections identified
        # by a keyword string. The header data are used by one or
        # more sections. The `section_header_map` maps each
        # section keyword to a specific header key.

        self._section_header_map = {}

        atoms_sections = ['Atoms', 'Velocities', 'Molecules']
        self._section_header_map.update(
            dict.fromkeys(atoms_sections, 'atoms'))

        bonds_sections = ['Bonds']
        self._section_header_map.update(
            dict.fromkeys(bonds_sections, 'bonds'))

        lines_sections = ['Lines']
        self._section_header_map.update(
            dict.fromkeys(lines_sections, 'lines'))

        ellipsoids_sections = ['Ellipsoids']
        self._section_header_map.update(
            dict.fromkeys(ellipsoids_sections, 'ellipsoids'))

        triangles_sections = ['Triangles']
        self._section_header_map.update(
            dict.fromkeys(triangles_sections, 'triangles'))

        bodies_sections = ['Bodies']
        self._section_header_map.update(
            dict.fromkeys(bodies_sections, 'bodies'))

        angles_sections = ['Angles']
        self._section_header_map.update(
            dict.fromkeys(angles_sections, 'angles'))

        dihedrals_sections = ['Dihedrals']
        self._section_header_map.update(
            dict.fromkeys(dihedrals_sections, 'dihedrals'))

        impropers_sections = ['Impropers']
        self._section_header_map.update(
            dict.fromkeys(impropers_sections, 'impropers'))

        atom_types_sections = ['Masses', 'Pair Coeffs']
        self._section_header_map.update(
            dict.fromkeys(atom_types_sections, 'atom types'))

        bond_types_sections = ['Bond Coeffs']
        self._section_header_map.update(
            dict.fromkeys(bond_types_sections, 'bond types'))

        angle_types_sections = \
            ['Angle Coeffs', 'BondBond Coeffs', 'BondAngle Coeffs']
        self._section_header_map.update(
            dict.fromkeys(angle_types_sections, 'angle types'))

        improper_types_sections = \
            ['AngleAngle Coeffs', 'Improper Coeffs']
        self._section_header_map.update(
            dict.fromkeys(improper_types_sections, 'improper types'))

        dihedral_types_sections = \
            ['Dihedral Coeffs', 'MiddleBondTorsion Coeffs',
             'EndBondTorsion Coeffs', 'AngleTorsion Coeffs',
             'AngleAngleTorsion Coeffs', 'BondBond13 Coeffs']
        self._section_header_map.update(
            dict.fromkeys(dihedral_types_sections, 'dihedral types'))

        self._properties['sections'] = self._section_header_map

    @property
    def atom_style(self):
        """Atom style."""
        return self._atom_style

    @property
    def header_keys(self):
        """List of header key names"""
        return self._header_keys

    @property
    def section_properties(self):
        """List of section properties"""
        return self._section_properties

    @property
    def section_syntax_dict(self):
        """Section syntax dictionary."""
        return self._section_syntax_dict
