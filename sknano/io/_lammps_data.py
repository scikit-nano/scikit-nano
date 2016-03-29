# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS data format (:mod:`sknano.io._lammps_data`)
====================================================================

.. currentmodule:: sknano.io._lammps_data

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
import os

import numpy as np

from monty.io import zopen
from sknano.core import get_fpath, minmax
from sknano.core.atoms import Atoms, MDAtoms, MDAtom as Atom
# from sknano.core.crystallography import Crystal3DLattice
from sknano.core.geometric_regions import Domain
from ._base import StructureData, StructureDataError, StructureDataFormatter, \
    StructureDataConverter, default_comment_line

__all__ = ['DATA', 'DATAData', 'DATAReader', 'DATAWriter', 'DATAFormatter',
           'DATAError', 'DATAIO', 'DATAIOReader', 'DATAIOWriter',
           'DATAIOError', 'DATAIOFormatter', 'DATAFormatSpec',
           'DATAConverter', 'DATA2XYZConverter', 'LAMMPSDATA',
           'LAMMPSDATAReader', 'LAMMPSDATAWriter', 'LAMMPSDATAFormatter',
           'LAMMPSDATAFormatSpec', 'LAMMPSDATAIOError',
           'LAMMPSDATA2XYZConverter', 'atom_styles', 'lammps_atom_styles']


class DATAReader(StructureData):
    """`StructureData` class for reading `LAMMPS data` file format.

    Parameters
    ----------
    fpath : str
        `LAMMPS data` file path
    atom_style : {'full', 'atomic'}, optional
    bond_style : {None, 'class2', 'fene', 'fene/expand', 'harmonic', 'morse', \
                  'nonlinear', 'quartic'}
    angle_style
    dihedral_style
    improper_style
    pair_style

    Attributes
    ----------
    domain : :class:`Domain`

    """
    def __init__(self, fpath, atom_style='full', bond_style=None,
                 angle_style=None, dihedral_style=None, improper_style=None,
                 pair_style=None, formatter=None, **kwargs):
        super().__init__(fpath=fpath, **kwargs)

        self.header_data = OrderedDict()
        self.section_data = OrderedDict()
        self.domain = Domain()

        if formatter is None or not isinstance(formatter, DATAFormatter):
            formatter = DATAFormatter(atom_style=atom_style,
                                      bond_style=bond_style,
                                      angle_style=angle_style,
                                      dihedral_style=dihedral_style,
                                      improper_style=improper_style,
                                      pair_style=pair_style)

        self.formatter = formatter
        self.section_attrs = formatter.section_attrs
        self.section_attrs_specs = formatter.section_attrs_specs
        self.fmtstr = "{fpath!r}, " + formatter.fmtstr

        if self.fpath is not None:
            self.read()

    def __str__(self):
        strrep = super().__str__()
        formatter = self.formatter
        if formatter is not None:
            strrep = '\n'.join((strrep, str(formatter)))
        return strrep

    @property
    def headers(self):
        """Alias for :attr:`~DATAReader.header_data`.

        Returns
        -------
        :class:`python:dict`
            :class:`python:dict` of dump file header values

        """
        return self.header_data

    @property
    def sections(self):
        """Alias for :attr:`~DATAReader.section_data`.

        Returns
        -------
        :class:`python:dict`
            :class:`python:dict` of dump file section data

        """
        return self.section_data

    def read(self):
        """Read data file."""
        self.structure.clear()
        try:
            with zopen(self.fpath) as f:
                self.comment_line = f.readline().strip()

                while True:
                    line = f.readline().strip()
                    if len(line) == 0:
                        continue
                    found = False
                    for header in list(header_specs.keys()):
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
                    for section, header in list(section_header_map.items()):
                        if section in line:
                            found = True
                            f.readline()
                            n = 0
                            data = []
                            while n < self.header_data[header]:
                                tmp = []
                                line = f.readline().strip().split()
                                for i, attrs in enumerate(
                                    self.section_attrs_specs[
                                        section].values()):
                                    try:
                                        tmp.append(
                                            attrs['dtype'](float(line[i])))
                                    except IndexError:
                                        break
                                data.append(tmp)
                                n += 1
                            self.section_data[section] = data[:]
                    f.readline()
                    line = f.readline().strip()
                    if len(line) == 0:
                        break
            self._parse_domain()
            self._parse_atoms()
            self._parse_atom_types()
            # self._parse_bonds()
            # self._parse_dihedrals()
            # self._parse_impropers()
            # self._parse_ellipsoids()
            # self._parse_lines()
            # self._parse_triangles()
            # self._parse_bodies()
        except (IOError, OSError) as e:
            print(e)

    def _parse_atoms(self):
        """Populate `Atoms` object with `Atom` objects"""
        try:
            atoms_section = self.section_data['Atoms']
            atoms_section_attrs = self.section_attrs['Atoms']
        except KeyError:
            atoms_section = []
            atoms_section_attrs = []

        try:
            masses_section = self.section_data['Masses']
            masses_section_attrs = self.section_attrs['Masses']
        except KeyError:
            masses_section = []
            masses_section_attrs = []

        try:
            velocities_section = self.section_data['Velocities']
            velocities_section_attrs = self.section_attrs['Velocities']
        except KeyError:
            velocities_section = []
            velocities_section_attrs = []

        for line in atoms_section:
            atom_kwargs = {}

            for kw in dir(Atom()):
                if kw in atoms_section_attrs:
                    try:
                        atom_kwargs[kw] = \
                            line[
                                self.section_attrs_specs['Atoms'][kw]['index']]
                    except IndexError:
                        pass
                elif kw in masses_section_attrs:
                    type = \
                        line[self.section_attrs_specs[
                            'Atoms']['type']['index']]
                    atom_kwargs[kw] = \
                        masses_section[type - 1][
                            self.section_attrs_specs['Masses'][kw]['index']]
                elif kw in velocities_section_attrs and \
                        len(velocities_section) == len(atoms_section):
                    id = line[self.section_attrs_specs['Atoms']['id']['index']]
                    for velocity in velocities_section:
                        atom_id = \
                            velocity[self.section_attrs_specs[
                                'Velocities']['id']['index']]
                        if atom_id == id:
                            atom_kwargs[kw] = \
                                velocity[self.section_attrs_specs[
                                    'Velocities'][kw]['index']]
                # else:
                #     print('unknown atom keyword: {}'.format(kw))

            atom = Atom(**atom_kwargs)
            self.atoms.append(atom)

    def _parse_atom_types(self):
        Ntypes = self.atoms.Ntypes
        typemap = self.atoms.typemap
        if Ntypes != self.header_data['atom types']:
            for atomtype in range(1, self.header_data['atom types'] + 1):
                if atomtype not in typemap:
                    try:
                        mass = self.section_data['Masses'][atomtype - 1][
                            self.section_attrs_specs['Masses']['mass']['index']
                        ]
                        self.atoms.add_type(Atom(type=atomtype, mass=mass))
                    except KeyError:
                        self.atoms.add_type(Atom(type=atomtype))

    def _parse_domain(self):
        domain = self.domain
        bounding_box = domain.bounding_box
        for dim in ('x', 'y', 'z'):
            bounds = \
                self.header_data[' '.join([dim + lim for lim in ('lo', 'hi')])]
            [setattr(bounding_box, dim + lim, value) for
             lim, value in zip(('min', 'max'), bounds)]
        domain.bounding_box = bounding_box

        tilt_factors = 'xy xz yz'
        if tilt_factors in self.headers:
            domain.triclinic = True
            [setattr(domain, tilt_factor, value) for tilt_factor, value
             in zip(tilt_factors.split(), self.headers[tilt_factors])]

    def _update_headers(self, from_atoms=None, from_bonds=None,
                        from_angles=None, from_dihedrals=None,
                        from_impropers=None):
        headers = OrderedDict()
        if from_atoms is not None:
            atoms = from_atoms
            headers['atoms'] = atoms.Natoms
            headers['atom types'] = atoms.Ntypes

        if from_bonds is not None:
            bonds = from_bonds
            headers['bonds'] = bonds.Nbonds
            headers['bond types'] = bonds.Ntypes

        if from_angles is not None:
            angles = from_angles
            headers['angles'] = angles.Nangles
            headers['angle types'] = angles.Ntypes

        if from_dihedrals is not None:
            dihedrals = from_dihedrals
            headers['dihedrals'] = dihedrals.Ndihedrals
            headers['dihedral types'] = dihedrals.Ntypes

        if from_impropers is not None:
            impropers = from_impropers
            headers['impropers'] = impropers.Nimpropers
            headers['improper types'] = impropers.Ntypes

        # TODO: adding bounding_box limits to headers
        self.header_data = headers

    def _update_sections(self, from_atoms=None):
        sections = OrderedDict()
        if from_atoms is not None:
            atoms = from_atoms
            typemap = OrderedDict(sorted(atoms.typemap.items()))
            sections['Masses'] = [[type, attrmap['mass']] for type, attrmap in
                                  typemap.items()]
            atoms_section_attrs = self.section_attrs['Atoms']
            sections['Atoms'] = \
                [[getattr(atom, attr) for attr in atoms_section_attrs]
                 for atom in atoms]
            velocities_section_attrs = self.section_attrs['Velocities']
            sections['Velocities'] = \
                [[getattr(atom, attr) for attr in velocities_section_attrs]
                 for atom in atoms]
        self.section_data = sections

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

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(self.formatter.todict())
        return attr_dict

LAMMPSDATAReader = DATAIOReader = DATAReader


class DATAWriter:
    """`StructureWriter` class for writing `LAMMPS data` file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, structure=None,
              atoms=None, bounding_box=None, allow_triclinic_box=False,
              atom_style='full', bond_style=None, angle_style=None,
              dihedral_style=None, improper_style=None, pair_style=None,
              **kwargs):
        """Write structure data to file.

        Parameters
        ----------
        fname : :class:`~python:str`, optional
            Output file name.
        outpath : :class:`~python:str`, optional
            Output file path.
        fpath : :class:`~python:str`, optional
            Full path (directory path + file name) to output data file.
        structure : :class:`~sknano.core.structures.BaseStructure`, optional
        atoms : :class:`~sknano.core.atoms.Atoms`, optional
            An :class:`~sknano.core.atoms.Atoms` instance.
        bounding_box : :class:`~python:dict`, optional
            If `None`, determined automatically from the `atoms` coordinates.
        allow_triclinic_box : :class:`~python:bool`, optional
        verbose : :class:`~python:bool`, optional
            verbose output

        """
        if structure is None and atoms is None:
            raise ValueError('Expected either `structure` or `atoms` object.')

        if structure is not None and atoms is None:
            atoms = structure.atoms

        if fpath is None:
            if 'datafile' in kwargs and fname is None:
                fname = kwargs.pop('datafile')
            fpath = get_fpath(fname=fname, ext='data', outpath=outpath,
                              overwrite=True, add_fnum=False)

        if not isinstance(atoms, MDAtoms):
            atoms = MDAtoms(atoms)

        formatter = DATAFormatter(atom_style=atom_style,
                                  bond_style=bond_style,
                                  angle_style=angle_style,
                                  dihedral_style=dihedral_style,
                                  improper_style=improper_style,
                                  pair_style=pair_style)

        data = DATAData(formatter=formatter)
        domain = data.domain
        if bounding_box is not None:
            domain.update(from_region=bounding_box,
                          allow_triclinic_box=allow_triclinic_box,
                          **kwargs)
        else:
            lattice = None
            if structure is not None and structure.lattice is not None:
                lattice = structure.lattice
            elif atoms.lattice is not None:
                lattice = atoms.lattice

            if lattice is not None:
                domain.update(from_lattice=lattice,
                              allow_triclinic_box=allow_triclinic_box,
                              **kwargs)
            else:
                domain.update(from_array=atoms.coords,
                              allow_triclinic_box=allow_triclinic_box,
                              **kwargs)

        data._update_headers(from_atoms=atoms)
        data._update_sections(from_atoms=atoms)
        data.write(datafile=fpath, atoms=atoms, **kwargs)

LAMMPSDATAWriter = DATAIOWriter = DATAWriter


class DATAData(DATAReader):
    """Class for reading and writing `StructureData` in `LAMMPS data` format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None, **kwargs):
        super().__init__(fpath, **kwargs)

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
                raise StructureDataError(e)
        attr_name = self.section_attrs[section][colidx]
        attr_dtype = self.section_attrs_specs[section][attr_name]['dtype']
        new_data = np.asarray(new_data, dtype=attr_dtype)

        for i, atom in enumerate(self.atoms):
            self.section_data[section][i][colidx] = \
                attr_dtype(float(new_data[i]))
            setattr(atom, attr_name, attr_dtype(float(new_data[i])))

    def viz(self, isnap):
        pass

    def write(self, datafile=None, atoms=None, comment_line=None, **kwargs):
        """Write data file.

        Parameters
        ----------
        datafile : :class:`~python:str`, optional
        atoms : :class:`~sknano.core.atoms.Atoms`, optional
        comment_line : :class:`~python:str`, optional
            A string written to the first line of `data` file. If `None`,
            then it is the current version string of scikit-nano.

        """
        try:
            kwargs.update(self.kwargs)

            if not datafile:
                if self.fpath is None:
                    error_msg = 'Invalid `datafile` {}'.format(datafile)
                    raise ValueError(error_msg)
                else:
                    datafile = self.fpath
            elif self.fpath is None:
                self.fpath = datafile

            if comment_line is None:
                comment_line = default_comment_line

            if atoms is not None and isinstance(atoms, Atoms):
                if not isinstance(atoms, MDAtoms):
                    atoms = MDAtoms(atoms)
                self._atoms = atoms

            super()._update_atoms(**kwargs)
            atoms = self._atoms
            atoms.assign_unique_ids()
            atoms.assign_unique_types()
            self._update_attr_fmtstr_widths()

            try:
                with zopen(datafile, 'wt') as stream:
                    self._write_header(stream, comment_line)
                    self._write_domain(stream)
                    [getattr(self, '_write_' + section.lower())(stream)
                     for section in self.sections.keys()]
            except OSError as e:
                print(e)

            self._atoms = self._atoms_copy

        except (TypeError, ValueError) as e:
            print(e)

    def _update_attr_fmtstr_widths(self):
        atoms = self.atoms
        attr_fmtstr_width['id'] = len(str(atoms.Natoms)) + 1
        attr_fmtstr_width['type'] = len(str(atoms.Ntypes)) + 1
        attr_fmtstr_width['mol'] = len(str(np.max(atoms.mols))) + 1
        attr_fmtstr_width['q'] = \
            len('{:f}'.format(np.max(np.abs(minmax(atoms.charges))))) + 2
        attr_fmtstr_width['mass'] = \
            len('{:f}'.format(np.max(atoms.masses))) + 4
        attr_fmtstr_width['x'] = \
            len('{:f}'.format(np.max(np.abs(minmax(atoms.x))))) + 2
        attr_fmtstr_width['y'] = \
            len('{:f}'.format(np.max(np.abs(minmax(atoms.y))))) + 2
        attr_fmtstr_width['z'] = \
            len('{:f}'.format(np.max(np.abs(minmax(atoms.z))))) + 2
        attr_fmtstr_width['ix'] = \
            len(str(np.max(np.abs(minmax(atoms.ix))))) + 2
        attr_fmtstr_width['iy'] = \
            len(str(np.max(np.abs(minmax(atoms.iy))))) + 2
        attr_fmtstr_width['iz'] = \
            len(str(np.max(np.abs(minmax(atoms.iz))))) + 2
        attr_fmtstr_width['vx'] = \
            len('{:f}'.format(np.max(np.abs(minmax(atoms.vx))))) + 2
        attr_fmtstr_width['vy'] = \
            len('{:f}'.format(np.max(np.abs(minmax(atoms.vy))))) + 2
        attr_fmtstr_width['vz'] = \
            len('{:f}'.format(np.max(np.abs(minmax(atoms.vz))))) + 2

        # attr_fmtstr_width['lx'] = \
        #     len('{:f}'.format(np.max(atoms.lx))) + 2
        # attr_fmtstr_width['ly'] = \
        #     len('{:f}'.format(np.max(atoms.ly))) + 2
        # attr_fmtstr_width['lz'] = \
        #     len('{:f}'.format(np.max(atoms.lz))) + 2
        # attr_fmtstr_width['wx'] = \
        #     len('{:f}'.format(np.max(atoms.wx))) + 2
        # attr_fmtstr_width['wy'] = \
        #     len('{:f}'.format(np.max(atoms.wy))) + 2
        # attr_fmtstr_width['wz'] = \
        #     len('{:f}'.format(np.max(atoms.wz))) + 2

        # attr_fmtstr_width['ervel'] = \
        #     len('{:f}'.format(np.max(atoms.ervel))) + 2
        # attr_fmtstr_width['shapex'] = \
        #     len('{:f}'.format(np.max(atoms.shapex))) + 2
        # attr_fmtstr_width['shapey'] = \
        #     len('{:f}'.format(np.max(atoms.shapey))) + 2
        # attr_fmtstr_width['shapez'] = \
        #     len('{:f}'.format(np.max(atoms.shapez))) + 2

        # attr_fmtstr_width['quatw'] = \
        #     len('{:f}'.format(np.max(atoms.quatw))) + 2
        # attr_fmtstr_width['quati'] = \
        #     len('{:f}'.format(np.max(atoms.quati))) + 2
        # attr_fmtstr_width['quatj'] = \
        #     len('{:f}'.format(np.max(atoms.quatj))) + 2
        # attr_fmtstr_width['quatk'] = \
        #     len('{:f}'.format(np.max(atoms.quatk))) + 2

        attr_fmtstr_width['atom1'] = attr_fmtstr_width['atom2'] = \
            attr_fmtstr_width['atom3'] = attr_fmtstr_width['atom4'] = \
            attr_fmtstr_width['id']

        for attr_specs in self.section_attrs_specs.values():
            for attr, specs in attr_specs.items():
                specs['width'] = attr_fmtstr_width[attr]

    def _write_header(self, stream, comment_line):
        stream.write('# {}\n\n'.format(comment_line.lstrip('#').strip()))
        for header, value in self.headers.items():
            if header in list(header_specs.keys())[-4:]:
                continue
            try:
                s = ' '.join(map(str, value[:] + list((header,))))
            except TypeError:
                s = ' '.join(map(str, list((value, header))))
            finally:
                stream.write('{}\n'.format(s))
        stream.write('\n')

    def _write_domain(self, stream):
        domain = self.domain
        bounding_box = domain.bounding_box
        lohi_width = 0
        lohi_fmtstr = '{:.10f} {:.10f}'
        for dim in ('x', 'y', 'z'):
            lohi_width = \
                max(lohi_width, len(lohi_fmtstr.format(
                    getattr(bounding_box, dim + 'min'),
                    getattr(bounding_box, dim + 'max'))) + 4)

        for dim in ('x', 'y', 'z'):
            stream.write('{}{dim}lo {dim}hi\n'.format(
                lohi_fmtstr.format(
                    getattr(bounding_box, dim + 'min'),
                    getattr(bounding_box, dim + 'max')).ljust(lohi_width),
                dim=dim))

        if domain.triclinic:
            stream.write('{xy:.10f} {xz:.10f} {yz:.10f} xy xz yz\n'.format(
                         xy=domain.xy, xz=domain.xz, yz=domain.yz))

    def _write_masses(self, stream):
        type_width = self.section_attrs_specs['Masses']['type']['width']
        stream.write('\nMasses\n\n')
        for type, mass in self.sections['Masses']:
            stream.write('{}{:.4f}\n'.format(
                         '{:d}'.format(type).ljust(type_width), mass))

    def _write_atoms(self, stream):
        section_attrs_specs_items = self.section_attrs_specs['Atoms'].items()
        stream.write('\nAtoms # {}\n\n'.format(self.formatter.atom_style))
        for atom in self.atoms:
            line = ''
            for attr, specs in section_attrs_specs_items:
                line += "{:>{}}".format(specs['fmtstr'].format(
                                        getattr(atom, attr)), specs['width'])
            line += '\n'
            stream.write(line)

    def _write_velocities(self, stream):
        section_attrs_specs_items = \
            self.section_attrs_specs['Velocities'].items()
        stream.write('\nVelocities\n\n')
        for atom in self.atoms:
            line = ''
            for attr, specs in section_attrs_specs_items:
                line += "{:>{}}".format(specs['fmtstr'].format(
                                        getattr(atom, attr)), specs['width'])
            line += '\n'
            stream.write(line)

    def _write_force_fields(self, stream):
        pass

    def _write_bonds(self, stream):
        pass

    def _write_angles(self, stream):
        pass

    def _write_dihedrals(self, stream):
        pass

    def _write_impropers(self, stream):
        pass

    @classmethod
    def formatter(cls, atom_style='full', **kwargs):
        """Return :class:`DATAFormatter` object."""
        return DATAFormatter(atom_style=atom_style, **kwargs)

DATA = LAMMPSDATA = DATAIO = DATAData


class DATAConverter(StructureDataConverter):
    """:class:`StructureDataConverter` class for converting `LAMMPS data`.

    Parameters
    ----------
    datafile : str

    """
    @property
    def datafile(self):
        """`LAMMPS data` file."""
        return self.infile


class DATA2XYZConverter(DATAConverter):
    """:class:`DATAConverter` class for converting to `xyz` format.

    .. versionadded:: 0.2.9

    Parameters
    ----------
    datafile : str

    """
    def __init__(self, datafile, **kwargs):
        xyzfile = os.path.splitext(datafile)[0] + '.xyz'
        super().__init__(infile=datafile, outfile=xyzfile, **kwargs)

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
        from ._xyz import XYZReader, XYZWriter

        kwargs.update(self.kwargs)

        datareader = DATAReader(self.infile, **kwargs)

        XYZWriter.write(fpath=self.outfile, atoms=datareader.atoms,
                        comment_line=datareader.comment_line, **kwargs)

        if return_reader:
            return XYZReader(self.outfile, **kwargs)

LAMMPSDATA2XYZConverter = DATA2XYZConverter


class DATAError(StructureDataError):
    """Exception class for :class:`DATAData` IO Errors."""
    pass

LAMMPSDATAIOError = DATAIOError = DATAError


class DATAFormatter(StructureDataFormatter):
    """`StructureDataFormatter` class the `LAMMPS data` format spec.

    Parameters
    ----------
    atom_style : str
        LAMMPS atom style.
    bond_style : str
    angle_style : str
    dihedral_style : str
    improper_style : str
    pair_style : str

    """
    def __init__(self, atom_style='full', bond_style=None, angle_style=None,
                 dihedral_style=None, improper_style=None, pair_style=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.atom_style = atom_style
        self.bond_style = bond_style
        self.angle_style = angle_style
        self.dihedral_style = dihedral_style
        self.improper_style = improper_style
        self.pair_style = pair_style

        self.section_attrs = OrderedDict()
        self.section_attrs['Masses'] = masses_section_attrs
        self.section_attrs['Atoms'] = atoms_section_attrs[atom_style]
        self.section_attrs['Velocities'] = velocities_section_attrs[atom_style]
        self.section_attrs['Bonds'] = bonds_section_attrs
        self.section_attrs['Angles'] = angles_section_attrs
        self.section_attrs['Dihedrals'] = dihedrals_section_attrs
        self.section_attrs['Impropers'] = impropers_section_attrs
        self.section_attrs['Pair Coeffs'] = \
            pair_coeffs_section_attrs[pair_style]
        self.section_attrs['Bond Coeffs'] = \
            bond_coeffs_section_attrs[bond_style]
        self.section_attrs['Angle Coeffs'] = \
            angle_coeffs_section_attrs[angle_style]
        self.section_attrs['Dihedral Coeffs'] = \
            dihedral_coeffs_section_attrs[dihedral_style]
        self.section_attrs['Improper Coeffs'] = \
            improper_coeffs_section_attrs[improper_style]
        # self.section_attrs['BondBond Coeffs'] = ['M', 'r1', 'r2']
        # self.section_attrs['BondAngle Coeffs'] = ['N1', 'N2', 'r1', 'r2']
        # self.section_attrs['MiddleBondTorsion Coeffs'] = \
        #     ['A1', 'A2', 'A3', 'r2']
        # self.section_attrs['EndBondTorsion Coeffs'] = \
        #     ['B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'r1', 'r3']
        # self.section_attrs['AngleTorsion Coeffs'] = \
        #     ['D1', 'D2', 'D3', 'E1', 'E2', 'E3', 'theta1', 'theta2']
        # self.section_attrs['AngleAngleTorsion Coeffs'] = \
        #     ['M', 'theta1', 'theta2']
        # self.section_attrs['BondBond13 Coeffs'] = ['N', 'r1', 'r3']
        # self.section_attrs['AngleAngle Coeffs'] = \
        #     ['M1', 'M2', 'M3', 'theta1', 'theta2', 'theta3']

        self.section_attrs_specs = OrderedDict()
        for section, attrs in list(self.section_attrs.items()):
            self.section_attrs_specs[section] = OrderedDict()
            for i, attr in enumerate(attrs):
                self.section_attrs_specs[section][attr] = \
                    {'dtype': attr_dtypes[attr] if attr in attr_dtypes
                     else float,
                     'colnum': i + 1,
                     'index': i,
                     'fmtstr': attr_fmtstr[attr] if attr in attr_fmtstr
                     else '{:f}',
                     'width': attr_fmtstr_width[attr] if attr in
                     attr_fmtstr_width else 14}

        self.fmtstr = "atom_style={atom_style!r}, " + \
            "bond_style={bond_style!r}, " + \
            "angle_style={angle_style!r}, " + \
            "dihedral_style={dihedral_style!r}, " + \
            "improper_style={improper_style!r}, " + \
            "pair_style={pair_style!r}"

    def __str__(self):
        strrep = super().__str__()
        items = ['atom_style', 'bond_style', 'dihedral_style',
                 'improper_style', 'pair_style']
        values = [self.atom_style, self.bond_style, self.dihedral_style,
                  self.improper_style, self.pair_style]
        table = list(zip(items, values))
        strrep = '\n'.join((strrep, table))
        return strrep

    @property
    def atom_style(self):
        """LAMMPS data atom_style."""
        return self._atom_style

    @atom_style.setter
    def atom_style(self, value):
        if value not in atom_styles:
            error_msg = "Unknown `atom_style`: {}\n".format(value) + \
                "Expected one of: {}".format(atom_styles)
            raise ValueError(error_msg)
        self._atom_style = value

    @property
    def angle_style(self):
        """LAMMPS data `angle_style`."""
        return self._angle_style

    @angle_style.setter
    def angle_style(self, value):
        if value not in angle_styles:
            error_msg = "Unknown `angle_style`: {}\n".format(value) + \
                "Expected one of: {}".format(angle_styles)
            raise ValueError(error_msg)
        self._angle_style = value

    @property
    def bond_style(self):
        """LAMMPS data bond_style."""
        return self._bond_style

    @bond_style.setter
    def bond_style(self, value):
        if value not in bond_styles:
            error_msg = "Unknown `bond_style`: {}\n".format(value) + \
                "Expected one of: {}".format(bond_styles)
            raise ValueError(error_msg)
        self._bond_style = value

    @property
    def dihedral_style(self):
        """LAMMPS data `dihedral_style`."""
        return self._dihedral_style

    @dihedral_style.setter
    def dihedral_style(self, value):
        if value not in dihedral_styles:
            error_msg = "Unknown `dihedral_style`: {}\n".format(value) + \
                "Expected one of: {}".format(dihedral_styles)
            raise ValueError(error_msg)
        self._dihedral_style = value

    @property
    def improper_style(self):
        """LAMMPS data `improper_style`."""
        return self._improper_style

    @improper_style.setter
    def improper_style(self, value):
        if value not in improper_styles:
            error_msg = "Unknown `improper_style`: {}\n".format(value) + \
                "Expected one of: {}".format(improper_styles)
            raise ValueError(error_msg)
        self._improper_style = value

    @property
    def pair_style(self):
        """LAMMPS data `pair_style`."""
        return self._pair_style

    @pair_style.setter
    def pair_style(self, value):
        if value not in pair_styles:
            error_msg = "Unknown `pair_style`: {}\n".format(value) + \
                "Expected one of: {}".format(pair_styles)
            raise ValueError(error_msg)
        self._pair_style = value

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(atom_style=self.atom_style,
                              angle_style=self.angle_style,
                              bond_style=self.bond_style,
                              dihedral_style=self.dihedral_style,
                              improper_style=self.improper_style,
                              pair_style=self.pair_style))
        return attr_dict

DATAFormatSpec = LAMMPSDATAFormatSpec = LAMMPSDATAFormatter = \
    DATAIOFormatter = DATAFormatter

header_specs = OrderedDict()
[header_specs.update({key: {'dtype': int, 'items': 1} for key in
                     ['atoms', 'atom types', 'bonds', 'bond types',
                      'angles', 'angle types', 'dihedrals', 'dihedral types',
                      'impropers', 'improper types',
                      'ellipsoids', 'lines', 'triangles', 'bodies',
                      'extra bond per atom', 'extra angle per atom',
                      'extra dihedral per atom', 'extra improper per atom',
                      'extra special per atom']})]
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

attr_dtypes = {'id': int, 'type': int, 'mol': int, 'q': float, 'mass': float,
               'x': float, 'y': float, 'z': float,
               'ix': int, 'iy': int, 'iz': int,
               'vx': float, 'vy': float, 'vz': float,
               'lx': float, 'ly': float, 'lz': float,
               'wx': float, 'wy': float, 'wz': float,
               'ervel': float,
               'shapex': float, 'shapey': float, 'shapez': float,
               'quatw': float, 'quati': float, 'quatj': float, 'quatk': float,
               'atom1': int, 'atom2': int, 'atom3': int, 'atom4': int}

attr_fmtstr = {key: '{:d}' if dtype == int else '{:f}'
               for key, dtype in attr_dtypes.items()}
attr_fmtstr_width = {key: 5 if dtype == int else 16
                     for key, dtype in attr_dtypes.items()}

masses_section_attrs = ['type', 'mass']

atoms_section_attrs = OrderedDict()
atoms_section_attrs['angle'] = ['id', 'mol', 'type', 'x', 'y', 'z']
atoms_section_attrs['atomic'] = ['id', 'type', 'x', 'y', 'z']
atoms_section_attrs['body'] = ['id', 'type', 'bodyflag', 'mass', 'x', 'y', 'z']
atoms_section_attrs['bond'] = ['id', 'mol', 'type', 'x', 'y', 'z']
atoms_section_attrs['charge'] = ['id', 'type', 'q', 'x', 'y', 'z']
atoms_section_attrs['dipole'] = \
    ['id', 'type', 'q', 'x', 'y', 'z', 'mux', 'muy', 'muz']
atoms_section_attrs['electron'] = \
    ['id', 'type', 'q', 'spin', 'eradius', 'x', 'y', 'z']
atoms_section_attrs['ellipsoid'] = \
    ['id', 'type', 'ellipsoidflag', 'density', 'x', 'y', 'z']
atoms_section_attrs['full'] = ['id', 'mol', 'type', 'q', 'x', 'y', 'z']
atoms_section_attrs['line'] = \
    ['id', 'mol', 'type', 'lineflag', 'density', 'x', 'y', 'z']
atoms_section_attrs['meso'] = ['id', 'type', 'rho', 'e', 'cv', 'x', 'y', 'z']
atoms_section_attrs['molecular'] = ['id', 'mol', 'type', 'x', 'y', 'z']
atoms_section_attrs['peri'] = \
    ['id', 'type', 'volume', 'density', 'x', 'y', 'z']
atoms_section_attrs['sphere'] = \
    ['id', 'type', 'diameter', 'density', 'x', 'y', 'z']
atoms_section_attrs['template'] = \
    ['id', 'mol', 'template-index', 'template-atom', 'type', 'x', 'y', 'z']
atoms_section_attrs['tri'] = \
    ['id', 'mol', 'type', 'triangleflag', 'density', 'x', 'y', 'z']
atoms_section_attrs['wavepacket'] = \
    ['id', 'type', 'charge', 'spin', 'eradius', 'etag', 'cs_re', 'cs_im',
     'x', 'y', 'z']
# atoms_section_attrs['hybrid'] = ['id', 'type', 'x', 'y', 'z', '...']

atom_styles = list(atoms_section_attrs.keys())
lammps_atom_styles = atom_styles

[attrs.extend(['ix', 'iy', 'iz']) for attrs in atoms_section_attrs.values()]

velocities_section_attrs = OrderedDict()
[velocities_section_attrs.update({atom_style: ['id', 'vx', 'vy', 'vz']})
 for atom_style in atom_styles]
velocities_section_attrs['electron'].append('ervel')
velocities_section_attrs['ellipsoid'].extend(['lx', 'ly', 'lz'])
velocities_section_attrs['sphere'].extend(['wx', 'wy', 'wz'])
# velocities_section_attrs['hybrid'].append('...')

bonds_section_attrs = ['id', 'type', 'atom1', 'atom2']
angles_section_attrs = ['id', 'type', 'atom1', 'atom2', 'atom3']
dihedrals_section_attrs = ['id', 'type', 'atom1', 'atom2', 'atom3', 'atom4']
impropers_section_attrs = ['id', 'type', 'atom1', 'atom2', 'atom3', 'atom4']
ellipsoids_section_attrs = ['id', 'shapex', 'shapey', 'shapez',
                            'quatw', 'quati', 'quatj', 'quatk']

pair_coeffs_section_attrs = OrderedDict()
pair_coeffs_section_attrs[None] = []
pair_styles = list(pair_coeffs_section_attrs.keys())

angle_coeffs_section_attrs = OrderedDict()
angle_coeffs_section_attrs[None] = []
angle_styles = list(angle_coeffs_section_attrs.keys())

bond_coeffs_section_attrs = OrderedDict()
bond_coeffs_section_attrs[None] = []
bond_coeffs_section_attrs['class2'] = ['type', 'R0', 'K2', 'K3', 'K4']
bond_coeffs_section_attrs['fene'] = \
    ['type', 'K', 'R0', 'epsilon', 'sigma']
bond_coeffs_section_attrs['fene/expand'] = \
    ['type', 'K', 'R0', 'epsilon', 'sigma', 'delta']
bond_coeffs_section_attrs['harmonic'] = ['type', 'K', 'r0']
bond_coeffs_section_attrs['morse'] = ['type', 'D', 'alpha', 'r0']
bond_coeffs_section_attrs['nonlinear'] = \
    ['type', 'epsilon', 'r0', 'lambda']
bond_coeffs_section_attrs['quartic'] = \
    ['type', 'K', 'B1', 'B2', 'Rc', 'U0']
bond_coeffs_section_attrs['table'] = ['type', 'filename', 'keyword']
bond_styles = list(bond_coeffs_section_attrs.keys())
bond_style_args = OrderedDict()
bond_style_args['table'] = ['style', 'N']

dihedral_coeffs_section_attrs = OrderedDict()
dihedral_coeffs_section_attrs[None] = []
dihedral_coeffs_section_attrs['charmm'] = \
    ['type', 'K', 'n', 'd', 'weight']
dihedral_coeffs_section_attrs['class2'] = []
dihedral_coeffs_section_attrs['harmonic'] = []
dihedral_coeffs_section_attrs['helix'] = []
dihedral_coeffs_section_attrs['multi/harmonic'] = []
dihedral_coeffs_section_attrs['opls'] = []
dihedral_styles = list(dihedral_coeffs_section_attrs.keys())

improper_coeffs_section_attrs = OrderedDict()
improper_coeffs_section_attrs[None] = []
improper_coeffs_section_attrs['class2'] = []
improper_coeffs_section_attrs['cvff'] = []
improper_coeffs_section_attrs['harmonic'] = []
improper_coeffs_section_attrs['umbrella'] = []
improper_styles = list(improper_coeffs_section_attrs.keys())
