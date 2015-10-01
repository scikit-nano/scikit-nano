# -*- coding: utf-8 -*-
"""
====================================================================
LAMMPS data format (:mod:`sknano.io._lammps_data_format`)
====================================================================

.. currentmodule:: sknano.io._lammps_data_format

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
import os

import numpy as np

from monty.io import zopen
from sknano.core import get_fpath
# from sknano.core.crystallography import Crystal3DLattice
from sknano.core.geometric_regions import generate_bounding_box, Cuboid
from ._base import Atom, StructureIO, StructureIOError, StructureConverter, \
    default_comment_line

__all__ = ['DATAReader', 'DATAWriter', 'DATAData', 'DATAFormatSpec',
           'DATAIOError', 'DATA2XYZConverter', 'LAMMPSDATAReader',
           'LAMMPSDATAWriter', 'LAMMPSDATA', 'LAMMPSDATAFormatSpec',
           'LAMMPSDATAIOError', 'LAMMPSDATA2XYZConverter',
           'atom_styles', 'lammps_atom_styles']


class Domain:
    pass


class DATAReader(StructureIO):
    """`StructureIO` class for reading `LAMMPS data` file format.

    Parameters
    ----------
    fpath : str
        `LAMMPS data` file path
    atom_style : {'full', 'atomic'}, optional

    """
    def __init__(self, fpath, atom_style='full', bond_style=None,
                 angle_style=None, dihedral_style=None, improper_style=None,
                 pair_style=None, **kwargs):
        super().__init__(fpath=fpath, **kwargs)

        self.header_data = OrderedDict()
        self.section_data = OrderedDict()
        self.domain = Domain()
        self.kwargs = OrderedDict()

        styles = dict(atom_style=atom_style, bond_style=bond_style,
                      angle_style=angle_style, dihedral_style=dihedral_style,
                      improper_style=improper_style, pair_style=pair_style)

        formatspec = DATAFormatSpec(**styles)
        self.atom_style = atom_style
        self.section_attrs = formatspec.section_attrs
        self.section_attrs_specs = formatspec.section_attrs_specs

        if self.fpath is not None:
            self.read()

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
        self.structure_data.clear()
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
            self._parse_bounding_box()
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

    def _parse_bounding_box(self):
        self.bounding_box = Cuboid()
        for dim in ('x', 'y', 'z'):
            bounds = \
                self.header_data[' '.join([dim + lim for lim in ('lo', 'hi')])]
            [setattr(self.bounding_box, dim + lim, value) for
             lim, value in zip(('min', 'max'), bounds)]

        self.kwargs['bounding_box'] = self.bounding_box

    def _parse_domain(self):
        tilt_factors = 'xy xz yz'
        if tilt_factors in self.headers:
            self.domain.triclinic = True
            [setattr(self.domain, tilt_factor, value) for tilt_factor, value
             in zip(tilt_factors.split(), self.headers[tilt_factors])]
        else:
            self.domain.triclinic = False

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


class DATAWriter:
    """`StructureWriter` class for writing `LAMMPS data` file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, structure=None,
              atoms=None, atom_style='full', bounding_box=None,
              comment_line=None, assert_unique_ids=False,
              enforce_consecutive_ids=True, pad_box=False,
              xpad=10., ypad=10., zpad=10., pad_tol=0.01,
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
        bounding_box : dict, optional
            If `None`, determined automatically from the `atoms` coordinates.
        comment_line : str, optional
            A string written to the first line of `data` file. If `None`,
            then it is set to the full path of the output `data` file.
        assert_unique_ids : bool, optional
            Check that each :class:`~sknano.core.atoms.Atom` in `atoms`
            has a unique :attr:`~sknano.core.atoms.Atom.id`.
            If the check fails, then assign a unique
            :attr:`~sknano.core.atoms.Atom.id`.
            to each :class:`~sknano.core.atoms.Atom`.
            If `assert_unique_ids` is True, but the id's are not unique,
            LAMMPS will not be able to read the data file.
        enforce_consecutive_ids : bool, optional
        pad_box : bool, optional
        xpad, ypad, zpad : float, optional
        pad_tol : float, optional
        verbose : bool, optional
            verbose output

        """
        if structure is None and atoms is None:
            raise ValueError('Expected either `structure` or `atoms` object.')

        if structure is not None and atoms is None:
            atoms = structure.atoms

        if fpath is None:
            fpath = get_fpath(fname=fname, ext='data', outpath=outpath,
                              overwrite=True, add_fnum=False)
        if comment_line is None:
            comment_line = default_comment_line

        atoms.rezero()
        atoms.assign_unique_types()
        typemap = atoms.typemap

        Natoms = atoms.Natoms
        Natoms_width = \
            8 if len(str(Natoms)) <= 12 else len(str(Natoms)) + 4
        Ntypes = atoms.Ntypes
        Ntypes_width = Natoms_width

        id_width = len(str(Natoms)) + 1
        type_width = len(str(Ntypes)) + 1

        if (enforce_consecutive_ids and
            atoms.ids.max() != atoms.Natoms) or \
                (not assert_unique_ids and
                 len(set(atoms.ids)) != atoms.Natoms):
            atoms.assign_unique_ids()

        if bounding_box is None:
            if structure is not None and structure.lattice is not None:
                bounding_box = \
                    generate_bounding_box(from_lattice=structure.lattice,
                                          center=atoms.centroid)
            else:
                bounding_box = \
                    generate_bounding_box(from_array=atoms.coords)

        if pad_box:
            boxpad = {'x': xpad, 'y': ypad, 'z': zpad}
            # for dim, pad in boxpad.items():
            for i, dim in enumerate(('x', 'y', 'z')):
                pad = boxpad[dim]
                dmin = dim + 'min'
                dmax = dim + 'max'
                if abs(getattr(bounding_box, dmin) -
                       atoms.coords[:, i].min()) < pad - pad_tol:
                    setattr(bounding_box, dmin,
                            getattr(bounding_box, dmin) - pad)
                if abs(getattr(bounding_box, dmax) -
                       atoms.coords[:, i].max()) < pad - pad_tol:
                    setattr(bounding_box, dmax,
                            getattr(bounding_box, dmax) + pad)

        if verbose:
            print('bounding_box: {}'.format(bounding_box))

        lohi_width = 0
        for dim in ('x', 'y', 'z'):
            lohi_width = \
                max(lohi_width, len('{:.6f} {:.6f}'.format(
                    getattr(bounding_box, dim + 'min'),
                    getattr(bounding_box, dim + 'max'))) + 4)

        with zopen(fpath, 'wt') as f:
            f.write('# {}\n\n'.format(comment_line.lstrip('#').strip()))

            f.write('{}atoms\n'.format(
                '{:d}'.format(Natoms).ljust(Natoms_width)))
            f.write('{}atom types\n\n'.format(
                '{:d}'.format(Ntypes).ljust(Ntypes_width)))

            for dim in ('x', 'y', 'z'):
                f.write('{}{dim}lo {dim}hi\n'.format(
                    '{:.6f} {:.6f}'.format(
                        getattr(bounding_box, dim + 'min'),
                        getattr(bounding_box, dim + 'max')).ljust(lohi_width),
                    dim=dim))

            f.write('\nMasses\n\n')
            for atomtype, properties in list(typemap.items()):
                f.write('{}{:.4f}\n'.format(
                    '{:d}'.format(atomtype).ljust(Natoms_width),
                    properties['mass']))

            f.write('\nAtoms\n\n')
            for atom in atoms:
                line = ''
                line += "{:>{}}".format(atom.id, id_width)
                line += "{:>{}}".format(atom.mol, 3)
                line += "{:>{}}".format(atom.type, type_width)
                line += "{:>{}}".format('{:.1f}'.format(atom.q), 4)
                line += "{:>{}}".format('{:f}'.format(atom.x), 14)
                line += "{:>{}}".format('{:f}'.format(atom.y), 14)
                line += "{:>{}}".format('{:f}'.format(atom.z), 14)
                line += "{:>{}}".format('{:d}'.format(atom.ix), 3)
                line += "{:>{}}".format('{:d}'.format(atom.iy), 3)
                line += "{:>{}}".format('{:d}'.format(atom.iz), 3)
                line += '\n'

                f.write(line)

            f.write('\nVelocities\n\n')
            for atom in atoms:
                line = ''
                line += "{:>{}}".format(atom.id, id_width)
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
            elif datafile is None or datafile == '':
                datafile = self.fpath

            # DATAWriter.write(fname=datafile, atoms=self.atoms,
            #                  comment_line=self.comment_line, **kwargs)

            self.atoms.rezero()
            self.atoms.assign_unique_types()
            self.atoms.assign_unique_ids()

            self._update_attr_fmtstr_widths()

            try:
                with zopen(datafile, 'wt') as fp:
                    self._write_header(fp)
                    self._write_bounding_box(fp)
                    [getattr(self, '_write_' + section.lower())(fp)
                     for section in self.sections.keys()]
            except OSError as e:
                print(e)

        except (TypeError, ValueError) as e:
            print(e)

    def _update_attr_fmtstr_widths(self):
        attr_fmtstr_width['id'] = len(str(self.atoms.Natoms)) + 1
        attr_fmtstr_width['type'] = len(str(self.atoms.Ntypes)) + 1
        attr_fmtstr_width['mol'] = len(str(np.max(self.atoms.mols))) + 1
        attr_fmtstr_width['q'] = \
            len('{:f}'.format(np.max(self.atoms.charges))) + 2
        attr_fmtstr_width['mass'] = \
            len('{:f}'.format(np.max(self.atoms.masses))) + 4
        attr_fmtstr_width['x'] = len('{:f}'.format(np.max(self.atoms.x))) + 2
        attr_fmtstr_width['y'] = len('{:f}'.format(np.max(self.atoms.y))) + 2
        attr_fmtstr_width['z'] = len('{:f}'.format(np.max(self.atoms.z))) + 2
        attr_fmtstr_width['ix'] = len(str(np.max(self.atoms.ix))) + 2
        attr_fmtstr_width['iy'] = len(str(np.max(self.atoms.iy))) + 2
        attr_fmtstr_width['iz'] = len(str(np.max(self.atoms.iz))) + 2
        attr_fmtstr_width['vx'] = len('{:f}'.format(np.max(self.atoms.vx))) + 2
        attr_fmtstr_width['vy'] = len('{:f}'.format(np.max(self.atoms.vy))) + 2
        attr_fmtstr_width['vz'] = len('{:f}'.format(np.max(self.atoms.vz))) + 2

        # attr_fmtstr_width['lx'] = \
        #     len('{:f}'.format(np.max(self.atoms.lx))) + 2
        # attr_fmtstr_width['ly'] = \
        #     len('{:f}'.format(np.max(self.atoms.ly))) + 2
        # attr_fmtstr_width['lz'] = \
        #     len('{:f}'.format(np.max(self.atoms.lz))) + 2
        # attr_fmtstr_width['wx'] = \
        #     len('{:f}'.format(np.max(self.atoms.wx))) + 2
        # attr_fmtstr_width['wy'] = \
        #     len('{:f}'.format(np.max(self.atoms.wy))) + 2
        # attr_fmtstr_width['wz'] = \
        #     len('{:f}'.format(np.max(self.atoms.wz))) + 2

        # attr_fmtstr_width['ervel'] = \
        #     len('{:f}'.format(np.max(self.atoms.ervel))) + 2
        # attr_fmtstr_width['shapex'] = \
        #     len('{:f}'.format(np.max(self.atoms.shapex))) + 2
        # attr_fmtstr_width['shapey'] = \
        #     len('{:f}'.format(np.max(self.atoms.shapey))) + 2
        # attr_fmtstr_width['shapez'] = \
        #     len('{:f}'.format(np.max(self.atoms.shapez))) + 2

        # attr_fmtstr_width['quatw'] = \
        #     len('{:f}'.format(np.max(self.atoms.quatw))) + 2
        # attr_fmtstr_width['quati'] = \
        #     len('{:f}'.format(np.max(self.atoms.quati))) + 2
        # attr_fmtstr_width['quatj'] = \
        #     len('{:f}'.format(np.max(self.atoms.quatj))) + 2
        # attr_fmtstr_width['quatk'] = \
        #     len('{:f}'.format(np.max(self.atoms.quatk))) + 2

        attr_fmtstr_width['atom1'] = attr_fmtstr_width['atom2'] = \
            attr_fmtstr_width['atom3'] = attr_fmtstr_width['atom4'] = \
            attr_fmtstr_width['id']

        for attr_specs in self.section_attrs_specs.values():
            for attr, specs in attr_specs.items():
                specs['width'] = attr_fmtstr_width[attr]

    def _write_header(self, fp):
        fp.write('# {}\n\n'.format(default_comment_line))
        for header, value in self.headers.items():
            if header in list(header_specs.keys())[-4:]:
                continue
            try:
                s = ' '.join(map(str, value[:] + list((header,))))
            except TypeError:
                s = ' '.join(map(str, list((value, header))))
            finally:
                fp.write('{}\n'.format(s))
        fp.write('\n')

    def _write_bounding_box(self, fp):
        lohi_width = 0
        lohi_fmtstr = '{:.10f} {:.10f}'
        for dim in ('x', 'y', 'z'):
            lohi_width = \
                max(lohi_width, len(lohi_fmtstr.format(
                    getattr(self.bounding_box, dim + 'min'),
                    getattr(self.bounding_box, dim + 'max'))) + 4)

        for dim in ('x', 'y', 'z'):
            fp.write('{}{dim}lo {dim}hi\n'.format(
                lohi_fmtstr.format(
                    getattr(self.bounding_box, dim + 'min'),
                    getattr(self.bounding_box, dim + 'max')).ljust(lohi_width),
                dim=dim))

        if self.domain.triclinic:
            fp.write('{xy:.10f} {xz:.10f} {yz:.10f} xy xz yz\n'.format(
                     xy=self.domain.xy, xz=self.domain.xz, yz=self.domain.yz))

    def _write_masses(self, fp):
        type_width = self.section_attrs_specs['Masses']['type']['width']
        fp.write('\nMasses\n\n')
        for type, mass in self.sections['Masses']:
            fp.write('{}{:.4f}\n'.format(
                     '{:d}'.format(type).ljust(type_width), mass))

    def _write_atoms(self, fp):
        fp.write('\nAtoms # {}\n\n'.format(self.atom_style))
        for atom in self.atoms:
            line = ''
            for attr, specs in self.section_attrs_specs['Atoms'].items():
                line += "{:>{}}".format(specs['fmtstr'].format(
                                        getattr(atom, attr)), specs['width'])
            line += '\n'
            fp.write(line)

    def _write_velocities(self, fp):
        fp.write('\nVelocities\n\n')
        for atom in self.atoms:
            line = ''
            for attr, specs in self.section_attrs_specs['Velocities'].items():
                line += "{:>{}}".format(specs['fmtstr'].format(
                                        getattr(atom, attr)), specs['width'])
            line += '\n'
            fp.write(line)

    def _write_force_fields(self, fp):
        pass

    def _write_bonds(self, fp):
        pass

    def _write_angles(self, fp):
        pass

    def _write_dihedrals(self, fp):
        pass

    def _write_impropers(self, fp):
        pass

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

        super().__init__(infile=self._datafile, outfile=self._xyzfile,
                         **kwargs)

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


class DATAFormatSpec:
    """`StructureFormatSpec` class the `LAMMPS data` format spec.

    Parameters
    ----------
    atom_style : {'full'}, optional
        LAMMPS atom style.

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
        self.section_attrs['Masses'] = ['type', 'mass']
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
        self.section_attrs['BondBond Coeffs'] = ['M', 'r1', 'r2']
        self.section_attrs['BondAngle Coeffs'] = ['N1', 'N2', 'r1', 'r2']
        self.section_attrs['MiddleBondTorsion Coeffs'] = \
            ['A1', 'A2', 'A3', 'r2']
        self.section_attrs['EndBondTorsion Coeffs'] = \
            ['B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'r1', 'r3']
        self.section_attrs['AngleTorsion Coeffs'] = \
            ['D1', 'D2', 'D3', 'E1', 'E2', 'E3', 'theta1', 'theta2']
        self.section_attrs['AngleAngleTorsion Coeffs'] = \
            ['M', 'theta1', 'theta2']
        self.section_attrs['BondBond13 Coeffs'] = ['N', 'r1', 'r3']
        self.section_attrs['AngleAngle Coeffs'] = \
            ['M1', 'M2', 'M3', 'theta1', 'theta2', 'theta3']

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

    @property
    def atom_style(self):
        return self._atom_style

    @atom_style.setter
    def atom_style(self, value):
        if value not in atom_styles:
            raise ValueError("Allowed `atom_style`'s:\n{}".format(
                list(atom_styles)))
        self._atom_style = value

    @atom_style.deleter
    def atom_style(self):
        del self._atom_style

LAMMPSDATAFormatSpec = DATAFormatSpec

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
               'atom1': int, 'atom2': int, 'atom3': int, 'atom4': int,
               'bond_id': int, 'bond_type': int,
               'angle_id': int, 'angle_type': int,
               'dihedral_id': int, 'dihedral_type': int,
               'improper_id': int, 'improper_type': int}

attr_fmtstr = {key: '{:d}' if dtype == int else '{:f}'
               for key, dtype in attr_dtypes.items()}
attr_fmtstr_width = {key: 5 if dtype == int else 16
                     for key, dtype in attr_dtypes.items()}

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

bonds_section_attrs = ['bond_id', 'bond_type', 'atom1', 'atom2']
angles_section_attrs = ['angle_id', 'angle_type', 'atom1', 'atom2', 'atom3']
dihedrals_section_attrs = ['dihedral_id', 'dihedral_type',
                           'atom1', 'atom2', 'atom3', 'atom4']
impropers_section_attrs = ['improper_id', 'improper_type',
                           'atom1', 'atom2', 'atom3', 'atom4']
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
bond_styles = list(bond_coeffs_section_attrs.keys())

dihedral_coeffs_section_attrs = OrderedDict()
dihedral_coeffs_section_attrs[None] = []
dihedral_styles = list(dihedral_coeffs_section_attrs.keys())

improper_coeffs_section_attrs = OrderedDict()
improper_coeffs_section_attrs[None] = []
improper_styles = list(improper_coeffs_section_attrs.keys())
