"""
===========================================================================
Structure readers (:mod:`sknano.structure_io._readers`)
===========================================================================

.. currentmodule:: sknano.structure_io._readers

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from abc import ABCMeta, abstractmethod

from pksci.chemistry import Atom, Atoms

__all__ = ['StructureReader', 'DATAReader', 'XYZReader']


class StructureReader(object):
    __metaclass__ = ABCMeta
    """Abstract superclass for reading structure data.

    Parameters
    ----------
    fname : str
        structure file

    """
    def __init__(self, fname=None):
        self._atoms = Atoms()
        self._fname = fname
        self._comment_line = None

    @property
    def atoms(self):
        return self._atoms

    @property
    def comment_line(self):
        return self._comment_line

    @property
    def fname(self):
        return self._fname

    @abstractmethod
    def _read(self):
        """Read in structure data from file"""
        return NotImplemented


class DATAReader(StructureReader):
    """Class for reading ``LAMMPS data`` file format.

    Parameters
    ----------
    datafile : str
        LAMMPS data file

    """
    def __init__(self, datafile, atom_style='full'):
        super(DATAReader, self).__init__(fname=datafile)

        from ._structure_specs import LAMMPSDATASpecs
        data_specs = LAMMPSDATASpecs(atom_style=atom_style)
        self._data_headers = data_specs.properties['headers']
        self._data_sections = data_specs.properties['sections']
        self._section_properties = data_specs.section_properties
        self._section_syntax_dict = data_specs.section_syntax_dict

        self._headers = {}
        self._sections = {}
        self._read()
        self._parse_atoms()
        self._Natoms = self._atoms.Natoms
        self._boxbounds = {}
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
                                enumerate(
                                    self._section_properties[
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
        """Populate Atoms list."""
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
        dict

        """
        try:
            section_data = self._sections[section_key]
            if colnum is not None:
                colidx = int(colnum - 1)
                return section_data[colidx]
            elif colname is not None:
                colidx = self._section_properties[section_key][colname]['index']
                return section_data[colidx]
            elif colindex is not None:
                colidx = int(colindex)
                return section_data[colindex]
        except (KeyError, TypeError, ValueError) as e:
            print(e)
        else:
            return section_data

    def map_colinfo(self):
        pass


class XYZReader(StructureReader):
    """Class for reading xyz chemical file format.

    Parameters
    ----------
    xyzfile : str
        xyz structure file

    """
    def __init__(self, xyzfile):
        super(XYZReader, self).__init__(fname=xyzfile)
        self._Natoms = None
        self._read()

    def _read(self):
        with open(self._fname, 'r') as f:
            self._Natoms = int(f.readline().strip())
            self._comment_line = f.readline().strip()
            lines = f.readlines()
            for line in lines:
                s = line.strip().split()
                if len(s) != 0:
                    atom = \
                        Atom(s[0], x=float(s[1]), y=float(s[2]), z=float(s[3]))
                    #atom.x, atom.y, atom.z = \
                    #    float(s[1]), float(s[2]), float(s[3])
                    self._atoms.append(atom)
