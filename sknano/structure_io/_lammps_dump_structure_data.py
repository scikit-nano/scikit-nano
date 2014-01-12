# -*- coding: utf-8 -*-
"""
==============================================================================
LAMMPS dump format (:mod:`sknano.structure_io._lammps_dump_structure_data`)
==============================================================================

.. currentmodule:: sknano.structure_io._lammps_dump_structure_data

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

import numpy as np

from pkshared.tools.fiofuncs import get_fpath

from ..chemistry import Atom, Atoms
from ._structure_data import StructureReader, StructureWriter, \
    StructureDataError


__all__ = ['DUMPReader', 'DUMPWriter', 'LAMMPSDUMP',
           'LAMMPSDUMPError']


class DUMPReader(StructureReader):
    """Class for reading ``LAMMPS dump`` file format.

    Parameters
    ----------
    fname : str
        LAMMPS dump file
    atom_style : {'full', 'atomic'}, optional

    """
    def __init__(self, fname=None, atom_style='full'):
        super(DUMPReader, self).__init__(fname=fname)

        from ._structure_specs import LAMMPSDUMPSpecs
        dump_specs = LAMMPSDUMPSpecs(atom_style=atom_style)
        self._dump_headers = dump_specs.properties['headers']
        self._dump_sections = dump_specs.properties['sections']
        self._section_properties = dump_specs.section_properties
        self._section_syntax_dict = dump_specs.section_syntax_dict

        self._headers = {}
        self._sections = {}
        self._boxbounds = {}

        if fname is not None:
            self.read()

    @property
    def headers(self):
        """DUMP file headers."""
        return self._headers

    @property
    def sections(self):
        """DUMP file sections."""
        return self._sections

    @property
    def boxbounds(self):
        """Box bounds."""
        return self._boxbounds

    def read(self):
        """Read dump file."""
        with open(self._fname, 'r') as f:
            self._comment_line = f.readline().strip()

            while True:
                line = f.readline().strip()
                if len(line) == 0:
                    continue
                found = False
                for key in self._dump_headers.iterkeys():
                    if key in line:
                        found = True
                        self._headers[key] = \
                            [self._dump_headers[key]['dtype'](s) for s in
                                [[ss for ss in line.split()][i] for i in
                                 range(self._dump_headers[key]['items'])]]
                        if len(self._headers[key]) == 1:
                            # since the list contains only one element,
                            # replace list value with element 0
                            self._headers[key] = self._headers[key][0]
                        break
                if not found:
                    break

            while True:
                found = False
                for section_key, header_key in self._dump_sections.iteritems():
                    if section_key in line:
                        found = True
                        f.readline()
                        Nitems = self._headers[header_key]
                        dumpdata = []
                        for n in xrange(Nitems):
                            tmp = []
                            line = f.readline().strip().split()
                            for i, props in \
                                enumerate(self._section_properties[
                                    section_key].itervalues()):
                                tmp.append(props['dtype'](line[i]))
                            dumpdata.append(tmp)
                        self._sections[section_key] = dumpdata[:]
                        break
                f.readline()
                line = f.readline().strip()
                if len(line) == 0:
                    break

        self._parse_atoms()
        self._parse_atomtypes()
        self._parse_boxbounds()

    def _parse_atoms(self):
        """Populate Atoms object with Atom objects"""
        atoms_section = self._sections['Atoms']
        atoms_section_syntax = self._section_syntax_dict['Atoms']

        masses_section = self._sections['Masses']
        masses_section_syntax = self._section_syntax_dict['Masses']

        velocities_section = self._sections['Velocities']
        velocities_section_syntax = self._section_syntax_dict['Velocities']

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
                else:
                    print('unknown atom keyword: {}'.format(kw))

            atom = Atom(**atom_kwargs)
            self._atoms.append(atom)

    def _parse_atomtypes(self):
        Ntypes = self._atoms.Ntypes
        atomtypes = self._atoms.atomtypes
        if Ntypes != self._headers['atom types']:
            for atomtype in xrange(1, self._headers['atom types'] + 1):
                if atomtype not in atomtypes:
                    mass = self._sections['Masses'][atomtype - 1][
                        self._section_properties['Masses']['mass']['index']]
                    self._atoms.add_atomtype(
                        Atom(atomtype=atomtype, mass=mass))

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
                except TypeError:
                    pass
        finally:
            return section_data


class DUMPWriter(StructureWriter):
    """Class for writing LAMMPS dump chemical file format."""

    @classmethod
    def write(cls, fname=None, atoms=None, boxbounds=None, comment_line=None,
              assume_unique_atoms=False, verbose=False):
        """Write structure dump to file.

        Parameters
        ----------
        fname : str
        atoms : `Atoms`
            An :py:class:`Atoms` instance.
        boxbounds : dict, optional
            If ``None``, determined automatically from atom coordinates.
        comment_line : str, optional
            A string written to the first line of ``dump`` file. If ``None``,
            then it is set to the full path of the output ``dump`` file.
        assume_unique_atoms : bool, optional
            Check that each Atom in Atoms has a unique atomID. If the check
            fails, then assign a unique atomID to each Atom.
            If ``assume_unique_atoms`` is True, but the atomID's are not
            unique, LAMMPS will not be able to read the dump file.
        verbose : bool, optional
            verbose output

        """
        if not isinstance(atoms, Atoms):
            raise TypeError('atoms argument must be an ``Atoms`` instance')
        else:
            fname = get_fpath(fname=fname, ext='dump', overwrite=True,
                              add_fnum=False)
            if comment_line is None:
                comment_line = fname

            atoms.fix_minus_zero_coords()

            atomtypes = atoms.atomtypes

            Natoms = atoms.Natoms
            Natoms_width = \
                8 if len(str(Natoms)) <= 12 else len(str(Natoms)) + 4
            Ntypes = atoms.Ntypes
            Ntypes_width = Natoms_width

            atomID_width = len(str(Natoms)) + 1
            atomtype_width = len(str(Ntypes)) + 1

            if not assume_unique_atoms and \
                    len(set(atoms.atom_ids)) != atoms.Natoms:
                for atomID, atom in enumerate(atoms, start=1):
                    atom.atomID = atomID

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
                    line += "{:>{}}\n".format('{:d}'.format(atom.nz), 3)

                    f.write(line)

                f.write('\nVelocities\n\n')
                for atom in atoms:
                    line = ''
                    line += "{:>{}}".format(atom.atomID, atomID_width)
                    line += "{:>{}}".format('{:f}'.format(atom.vx), 14)
                    line += "{:>{}}".format('{:f}'.format(atom.vy), 14)
                    line += "{:>{}}\n".format('{:f}'.format(atom.vz), 14)

                    f.write(line)


class aselect(object):
    pass


class tselect(object):
    pass


class LAMMPSDUMP(DUMPReader):
    """Class for reading and writing structure data in LAMMPS dump format.

    Parameters
    ----------
    fname : str, optional

    """
    def __init__(self, fname=None):
        super(LAMMPSDUMP, self).__init__(fname=fname)
        self._snaps = []
        self._nsnaps = self._nselect = 0
        self._names = {}
        self._aselect = aselect(self)
        self._tselect = tselect(self)

    def atom(self, n, l=[]):
        pass

    def read_all(self):
        pass

    def next(self):
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

        # for backwards compatibility with the pizza.py dump module,
        # first check positional arguments to see if this method was called
        # using the pizza.py dump module signature which expects positional
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
                raise LAMMPSDUMPError('replace called with invalid arguments')
        atom_attr = self._section_syntax_dict[section_key][colidx]
        attr_dtype = \
            self._section_properties[section_key][atom_attr]['dtype']
        new_data = np.asarray(new_data, dtype=attr_dtype)

        for i, atom in enumerate(self._atoms):
            self._sections[section_key][i][colidx] = attr_dtype(new_data[i])
            setattr(atom, atom_attr, attr_dtype(new_data[i]))

    def write(self, dumpfile=None):
        """Write dump file.

        Parameters
        ----------
        dumpfile : {None, str}, optional

        """
        try:
            if (dumpfile is None or dumpfile == '') and \
                    (self.fname is None or self.fname == ''):
                error_msg = '`dumpfile` must be a string at least 1 ' + \
                    'character long.'
                if dumpfile is None:
                    raise TypeError(error_msg)
                else:
                    raise ValueError(error_msg)
            else:
                dumpfile = self._fname
            DUMPWriter.write(fname=dumpfile, atoms=self._atoms,
                             boxbounds=self._boxbounds,
                             comment_line=self._comment_line)
        except (TypeError, ValueError) as e:
            print(e)


class LAMMPSDUMPError(StructureDataError):
    """Exception raised for failed method calls.

    Parameters
    ----------
    msg : str
        Error message

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)
