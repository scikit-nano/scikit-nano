# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS data format (:mod:`sknano.io._lammps_data_format`)
====================================================================

.. currentmodule:: sknano.io._lammps_data_format

"""
from __future__ import absolute_import, division, print_function
import six
from six.moves import range
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
import os

import numpy as np

from sknano.core import get_fpath
from ._base import Atom, StructureIO, StructureIOError, StructureConverter, \
    default_comment_line

__all__ = ['DATAReader', 'DATAWriter', 'DATAData', 'DATAFormatSpec',
           'DATAIOError', 'DATA2XYZConverter', 'LAMMPSDATAReader',
           'LAMMPSDATAWriter', 'LAMMPSDATA', 'LAMMPSDATAFormatSpec',
           'LAMMPSDATAIOError', 'LAMMPSDATA2XYZConverter',
           'atom_styles', 'lammps_atom_styles']


class DATAReader(StructureIO):
    """`StructureIO` class for reading `LAMMPS data` file format.

    Parameters
    ----------
    fpath : str
        `LAMMPS data` file path
    atom_style : {'full', 'atomic'}, optional

    """
    def __init__(self, fpath, atom_style='full', **kwargs):
        super(DATAReader, self).__init__(fpath=fpath, **kwargs)

        self.header_data = {}
        self.section_data = {}
        self.boxbounds = {}

        formatspec = DATAFormatSpec(atom_style=atom_style, **kwargs)
        self.section_attrs = formatspec.section_attrs
        self.section_attrs_specs = formatspec.section_attrs_specs

        if self.fpath is not None:
            self.read()

    @property
    def headers(self):
        """:class:`python:dict` of dump file header values"""
        return self.header_data

    @property
    def sections(self):
        return self.section_data

    def read(self):
        """Read data file."""
        self.structure_data.clear()
        try:
            with open(self.fpath, 'r') as f:
                self.comment_line = f.readline().strip()

                while True:
                    line = f.readline().strip()
                    if len(line) == 0:
                        continue
                    found = False
                    for header in six.iterkeys(header_specs):
                        if header in line:
                            found = True
                            self.header_data[header] = \
                                [header_specs[header]['dtype'](float(s))
                                 for s in [[ss for ss in line.split()][i]
                                           for i in range(header_specs[
                                               header]['items'])]]
                            if len(self.header_data[header]) == 1:
                                # if the list contains only one element,
                                # replace list with the first element
                                self.header_data[header] = \
                                    self.header_data[header][0]
                    if not found:
                        break

                while True:
                    found = False
                    for section, header in six.iteritems(section_header_map):
                        if section in line:
                            found = True
                            f.readline()
                            Nitems = self.header_data[header]
                            data = []
                            for n in range(Nitems):
                                tmp = []
                                line = f.readline().strip().split()
                                for i, attrs in enumerate(
                                    six.itervalues(self.section_attrs_specs[
                                        section])):
                                    try:
                                        tmp.append(
                                            attrs['dtype'](float(line[i])))
                                    except IndexError:
                                        break
                                data.append(tmp)
                            self.section_data[section] = data[:]
                    f.readline()
                    line = f.readline().strip()
                    if len(line) == 0:
                        break
            self._parse_atoms()
            self._parse_atomtypes()
            self._parse_boxbounds()
        except (IOError, OSError) as e:
            print(e)

    def _parse_atoms(self):
        """Populate `Atoms` object with `Atom` objects"""
        try:
            atoms_section = self.section_data['Atoms']
            atoms_section_attrs = self.section_attrs['Atoms']
        except KeyError:
            atoms_section = []
            atoms_section_attrs = {}

        try:
            masses_section = self.section_data['Masses']
            masses_section_attrs = self.section_attrs['Masses']
        except KeyError:
            masses_section = []
            masses_section_attrs = {}

        try:
            velocities_section = self.section_data['Velocities']
            velocities_section_attrs = self.section_attrs['Velocities']
        except KeyError:
            velocities_section = []
            velocities_section_attrs = {}

        atom_kwargs = {'atomID': None, 'moleculeID': None,
                       'q': None, 'atomtype': None, 'mass': None,
                       'x': None, 'y': None, 'z': None,
                       'vx': None, 'vy': None, 'vz': None}

        for line in atoms_section:
            for kw in six.iterkeys(atom_kwargs):
                if kw in atoms_section_attrs:
                    atom_kwargs[kw] = \
                        line[self.section_attrs_specs['Atoms'][kw]['index']]
                elif kw in masses_section_attrs:
                    atomtype = \
                        line[self.section_attrs_specs[
                            'Atoms']['atomtype']['index']]
                    atom_kwargs[kw] = \
                        masses_section[atomtype-1][
                            self.section_attrs_specs['Masses'][kw]['index']]
                elif kw in velocities_section_attrs and \
                        len(velocities_section) == len(atoms_section):
                    atomID = \
                        line[self.section_attrs_specs[
                            'Atoms']['atomID']['index']]
                    for velocity in velocities_section:
                        velocity_atomID = \
                            velocity[self.section_attrs_specs[
                                'Velocities']['atomID']['index']]
                        if velocity_atomID == atomID:
                            atom_kwargs[kw] = \
                                velocity[self.section_attrs_specs[
                                    'Velocities'][kw]['index']]
                #else:
                #    print('unknown atom keyword: {}'.format(kw))

            atom = Atom(**atom_kwargs)
            self.atoms.append(atom)

    def _parse_atomtypes(self):
        Ntypes = self.atoms.Ntypes
        atomtypes = self.atoms.atomtypes
        if Ntypes != self.header_data['atom types']:
            for atomtype in range(1, self.header_data['atom types'] + 1):
                if atomtype not in atomtypes:
                    mass = self.section_data['Masses'][atomtype - 1][
                        self.section_attrs_specs['Masses']['mass']['index']]
                    self.atoms.add_atomtype(
                        Atom(atomtype=atomtype, mass=mass))

    def _parse_boxbounds(self):
        for dim in ('x', 'y', 'z'):
            bounds = \
                self.header_data[' '.join([dim + lim for lim in ('lo', 'hi')])]
            self.boxbounds[dim] = {'min': bounds[0], 'max': bounds[-1]}

        self.kwargs['boxbounds'] = self.boxbounds

    def get(self, section, colnum=None, colname=None, colindex=None):
        """Return section with `section`.

        Parameters
        ----------
        section : str
        colnum : int, optional
        colname : str, optional
        colindex : int, optional

        Returns
        -------
        sequence

        """
        section_data = None
        try:
            section_data = self.section_data[section]
            section_attrs = self.section_attrs[section]
        except KeyError as e:
            print(e)
        else:
            try:
                colidx = None
                if colnum is not None:
                    colidx = int(colnum - 1)
                elif colname is not None:
                    colidx = \
                        self.section_attrs_specs[section][colname]['index']
                elif colindex is not None:
                    colidx = int(colindex)
            except (KeyError, TypeError, ValueError) as e:
                print(e)
            else:
                try:
                    colname = section_attrs[colidx]
                    coltype = \
                        self.section_attrs_specs[section][colname]['dtype']
                    section_data = \
                        np.asarray(
                            section_data, dtype=coltype)[:, colidx].tolist()
                except TypeError:
                    pass
        finally:
            return section_data

LAMMPSDATAReader = DATAReader


class DATAWriter(object):
    """`StructureWriter` class for writing `LAMMPS data` file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, atoms=None,
              atom_style='full', boxbounds=None, comment_line=None,
              assert_unique_ids=False, enforce_consecutive_atomIDs=True,
              pad_box=True, xpad=10., ypad=10., zpad=10., pad_tol=0.01,
              verbose=False, **kwargs):
        """Write structure data to file.

        Parameters
        ----------
        fname : str, optional
            Output file name.
        outpath : str, optional
            Output file path.
        fpath : str, optional
            Full path (directory path + file name) to output data file.
        atoms : :class:`~sknano.core.atoms.Atoms`
            An :class:`~sknano.core.atoms.Atoms` instance.
        boxbounds : dict, optional
            If `None`, determined automatically from the `atoms` coordinates.
        comment_line : str, optional
            A string written to the first line of `data` file. If `None`,
            then it is set to the full path of the output `data` file.
        assert_unique_ids : bool, optional
            Check that each :class:`~sknano.core.atoms.Atom` in `atoms`
            has a unique :attr:`~sknano.core.atoms.Atom.atomID`.
            If the check fails, then assign a unique
            :attr:`~sknano.core.atoms.Atom.atomID`.
            to each :class:`~sknano.core.atoms.Atom`.
            If `assert_unique_ids` is True, but the atomID's are not unique,
            LAMMPS will not be able to read the data file.
        enforce_consecutive_atomIDs : bool, optional
        pad_box : bool, optional
        xpad, ypad, zpad : float, optional
        pad_tol : float, optional
        verbose : bool, optional
            verbose output

        """
        if fpath is None:
            fpath = get_fpath(fname=fname, ext='data', outpath=outpath,
                              overwrite=True, add_fnum=False)
        if comment_line is None:
            comment_line = default_comment_line

        atoms.rezero()

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
                (not assert_unique_ids and
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
        if pad_box:
            #for dim, pad in boxpad.iteritems():
            for i, dim in enumerate(('x', 'y', 'z')):
                pad = boxpad[dim]
                if abs(boxbounds[dim]['min'] - atoms.coords[:, i].min()) \
                        < pad - pad_tol:
                    boxbounds[dim]['min'] = boxbounds[dim]['min'] - pad
                if abs(boxbounds[dim]['max'] - atoms.coords[:, i].max()) \
                        < pad - pad_tol:
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
            for atomtype, properties in six.iteritems(atomtypes):
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

LAMMPSDATAWriter = DATAWriter


class DATAData(DATAReader):
    """Class for reading and writing `StructureIO` in `LAMMPS data` format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None, **kwargs):
        super(DATAData, self).__init__(fpath, **kwargs)

    def delete(self, key):
        try:
            del self.header_data[key]
        except KeyError:
            try:
                del self.section_data[key]
            except KeyError as e:
                print(e)

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

    def replace(self, section, new_data, colnum=None,
                colname=None, colindex=None):
        """Replace section data.

        Parameters
        ----------
        section : str
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
                        self.section_attrs_specs[section][colname]['index']
                elif colindex is not None:
                    colidx = int(colindex)
            except (KeyError, TypeError, ValueError) as e:
                raise DATAIOError(e)
        attr_name = self.section_attrs[section][colidx]
        attr_dtype = self.section_attrs_specs[section][attr_name]['dtype']
        new_data = np.asarray(new_data, dtype=attr_dtype)

        for i, atom in enumerate(self.atoms):
            self.section_data[section][i][colidx] = \
                attr_dtype(float(new_data[i]))
            setattr(atom, attr_name, attr_dtype(float(new_data[i])))

    def viz(self, isnap):
        pass

    def write(self, datafile=None, **kwargs):
        """Write data file.

        Parameters
        ----------
        datafile : {None, str}, optional

        """
        try:
            kwargs.update(self.kwargs)

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

            DATAWriter.write(fname=datafile, atoms=self.atoms,
                             comment_line=self.comment_line, **kwargs)

        except (TypeError, ValueError) as e:
            print(e)

    @classmethod
    def format_spec(cls, atom_style='full', **kwargs):
        return DATAFormatSpec(atom_style=atom_style, **kwargs)

LAMMPSDATA = DATAData


class DATA2XYZConverter(StructureConverter):
    """
    `StructureConverter` class for converting `LAMMPS data` to `xyz` format.

    .. versionadded:: 0.2.9

    Parameters
    ----------
    datafile : str

    """
    def __init__(self, datafile, **kwargs):
        self._datafile = datafile
        self._xyzfile = os.path.splitext(self._datafile)[0] + '.xyz'

        super(DATA2XYZConverter, self).__init__(
            infile=self._datafile, outfile=self._xyzfile, **kwargs)

    @property
    def datafile(self):
        """`LAMMPS data` file."""
        return self.infile

    @property
    def xyzfile(self):
        """`xyz` file name."""
        return self.outfile

    def convert(self, return_reader=False, **kwargs):
        """Convert `LAMMPS data` to `xyz` chemical file format.

        Parameters
        ----------
        return_reader : bool, optional
            Return an instance of `XYZReader`

        Returns
        -------
        `XYZReader` (only if `return_reader` is True)

        """
        from ._xyz_format import XYZReader, XYZWriter

        kwargs.update(self.kwargs)

        datareader = DATAReader(self.infile, **kwargs)

        XYZWriter.write(fpath=self.outfile, atoms=datareader.atoms,
                        comment_line=datareader.comment_line, **kwargs)

        if return_reader:
            return XYZReader(self.outfile, **kwargs)

LAMMPSDATA2XYZConverter = DATA2XYZConverter


class DATAIOError(StructureIOError):
    pass

LAMMPSDATAIOError = DATAIOError


class DATAFormatSpec(object):
    """`StructureFormatSpec` class the `LAMMPS data` format spec.

    Parameters
    ----------
    atom_style : {'full'}, optional
        LAMMPS atom style.

    """

    def __init__(self, atom_style='full', bond_style=None, angle_style=None,
                 dihedral_style=None, improper_style=None, pair_style=None,
                 **kwargs):
        super(DATAFormatSpec, self).__init__(**kwargs)
        self.atom_style = atom_style
        self.bond_style = bond_style
        self.angle_style = angle_style
        self.dihedral_style = dihedral_style
        self.improper_style = improper_style
        self.pair_style = pair_style

        self.section_attrs = \
            {'Atoms': atoms_section_attrs[self.atom_style],
             'Masses': ['atomtype', 'mass'],
             'Velocities': velocities_section_attrs[self.atom_style],
             'Bonds': ['bond_id', 'bondtype', 'atom1', 'atom2'],
             'Dihedrals': ['dihedral_id', 'dihedraltype',
                           'atom1', 'atom2', 'atom3', 'atom4'],
             'Ellipsoids': ['atomID', 'shapex', 'shapey', 'shapez',
                            'quatw', 'quati', 'quatj', 'quatk']}

        self.section_attrs_specs = OrderedDict()
        for section, attrs in six.iteritems(self.section_attrs):
            self.section_attrs_specs[section] = OrderedDict()
            for i, attr in enumerate(attrs):
                self.section_attrs_specs[section][attr] = \
                    {'dtype': attr_dtypes[attr], 'colnum': i+1, 'index': i}

    @property
    def atom_style(self):
        return self._atom_style

    @atom_style.setter
    def atom_style(self, value):
        if value not in atom_styles:
            raise ValueError("Allowed `atom_style`'s:\n{}".format(
                list(atom_styles.keys())))
        self._atom_style = value

    @atom_style.deleter
    def atom_style(self):
        del self._atom_style

LAMMPSDATAFormatSpec = DATAFormatSpec

header_specs = OrderedDict()
header_specs.update(dict.fromkeys(['atoms', 'bonds', 'angles',
                                   'dihedrals', 'impropers',
                                   'atom types', 'bond types',
                                   'angle types', 'dihedral types',
                                   'improper types',
                                   'extra bond per atom',
                                   'extra angle per atom',
                                   'extra dihedral per atom',
                                   'extra improper per atom',
                                   'extra special per atom',
                                   'ellipsoids', 'lines',
                                   'triangles', 'bodies'],
                                  {'dtype': int, 'items': 1}))
header_specs.update(dict.fromkeys(['xlo xhi', 'ylo yhi', 'zlo zhi'],
                                  {'dtype': float, 'items': 2}))
header_specs.update(dict.fromkeys(['xy xz yz'],
                                  {'dtype': float, 'items': 3}))

# A LAMMPS data file is partitioned into sections identified
# by a keyword string. The header data are used by one or
# more sections. The `sections` dictionary maps each
# section keyword to a specific header key.

section_header_map = OrderedDict()
section_header_map.update(dict.fromkeys(['Atoms', 'Velocities', 'Molecules'],
                                        'atoms'))
section_header_map.update(dict.fromkeys(['Bonds'], 'bonds'))
section_header_map.update(dict.fromkeys(['Lines'], 'lines'))
section_header_map.update(dict.fromkeys(['Ellipsoids'], 'ellipsoids'))
section_header_map.update(dict.fromkeys(['Triangles'], 'triangles'))
section_header_map.update(dict.fromkeys(['Bodies'], 'bodies'))
section_header_map.update(dict.fromkeys(['Angles'], 'angles'))
section_header_map.update(dict.fromkeys(['Dihedrals'], 'dihedrals'))
section_header_map.update(dict.fromkeys(['Impropers'], 'impropers'))
section_header_map.update(dict.fromkeys(['Masses', 'Pair Coeffs'],
                                        'atom types'))
section_header_map.update(dict.fromkeys(['Bond Coeffs'], 'bond types'))
section_header_map.update(dict.fromkeys(['Angle Coeffs', 'BondBond Coeffs',
                                         'BondAngle Coeffs'], 'angle types'))
section_header_map.update(dict.fromkeys(['AngleAngle Coeffs',
                                         'Improper Coeffs'],
                                        'improper types'))
section_header_map.update(dict.fromkeys(['Dihedral Coeffs',
                                         'MiddleBondTorsion Coeffs',
                                         'EndBondTorsion Coeffs',
                                         'AngleTorsion Coeffs',
                                         'AngleAngleTorsion Coeffs',
                                         'BondBond13 Coeffs'],
                                        'dihedral types'))

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

atom_styles = {}
atom_styles['angle'] = ['atom-ID', 'molecule-ID', 'atom-type', 'x', 'y', 'z']
atom_styles['atomic'] = ['atom-ID', 'atom-type', 'x', 'y', 'z']
atom_styles['body'] = \
    ['atom-ID', 'atom-type', 'bodyflag', 'mass', 'x', 'y', 'z']
atom_styles['bond'] = ['atom-ID', 'molecule-ID', 'atom-type', 'x', 'y', 'z']
atom_styles['charge'] = ['atom-ID', 'atom-type', 'q', 'x', 'y', 'z']
atom_styles['dipole'] = \
    ['atom-ID', 'atom-type', 'q', 'x', 'y', 'z', 'mux', 'muy', 'muz']
atom_styles['electron'] = \
    ['atom-ID', 'atom-type', 'q', 'spin', 'eradius', 'x', 'y', 'z']
atom_styles['ellipsoid'] = \
    ['atom-ID', 'atom-type', 'ellipsoidflag', 'density', 'x', 'y', 'z']
atom_styles['full'] = \
    ['atom-ID', 'molecule-ID', 'atom-type', 'q', 'x', 'y', 'z']
atom_styles['line'] = \
    ['atom-ID', 'molecule-ID', 'atom-type', 'lineflag', 'density',
     'x', 'y', 'z']
atom_styles['meso'] = ['atom-ID', 'atom-type', 'rho', 'e', 'cv', 'x', 'y', 'z']
atom_styles['molecular'] = \
    ['atom-ID', 'molecule-ID', 'atom-type', 'x', 'y', 'z']
atom_styles['peri'] = \
    ['atom-ID', 'atom-type', 'volume', 'density', 'x', 'y', 'z']
atom_styles['sphere'] = \
    ['atom-ID', 'atom-type', 'diameter', 'density', 'x', 'y', 'z']
atom_styles['template'] = \
    ['atom-ID', 'molecule-ID', 'template-index', 'template-atom', 'atom-type',
     'x', 'y', 'z']
atom_styles['tri'] = \
    ['atom-ID', 'molecule-ID', 'atom-type', 'triangleflag', 'density',
     'x', 'y', 'z']
atom_styles['wavepacket'] = \
    ['atom-ID', 'atom-type', 'charge', 'spin', 'eradius', 'etag',
     'cs_re', 'cs_im', 'x', 'y', 'z']
#atom_styles['hybrid'] = ['atom-ID', 'atom-type', 'x', 'y', 'z', '...']

lammps_atom_styles = atom_styles

atoms_section_attrs = {}
for atom_style, attrs in six.iteritems(atom_styles):
    atom_style_attrs = atoms_section_attrs[atom_style] = []
    for attr in attrs:
        atom_style_attrs.append(attr.replace('-', ''))
    atom_style_attrs.extend(['nx', 'ny', 'nz'])

velocities_section_attrs = {}
velocities_section_attrs.update(dict.fromkeys(list(atom_styles.keys()),
                                              ['atomID', 'vx', 'vy', 'vz']))
velocities_section_attrs['electron'].append('ervel')
velocities_section_attrs['ellipsoid'].extend(['lx', 'ly', 'lz'])
velocities_section_attrs['sphere'].extend(['wx', 'wy', 'wz'])
#velocities_section_attrs['hybrid'].append('...')
