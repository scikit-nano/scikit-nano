# -*- coding: utf-8 -*-
"""
==============================================================================
Structure format specifications (:mod:`sknano.structure_io._structure_specs`)
==============================================================================

.. currentmodule:: sknano.structure_io._structure_specs

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from collections import OrderedDict

__all__ = ['StructureSpecs', 'LAMMPSDATASpecs']


class StructureSpecs(object):
    """Base class defining common properties for structure formats."""

    def __init__(self):
        self._properties = OrderedDict()

    @property
    def properties(self):
        """OrderedDict of format properties."""
        return self._properties


class LAMMPSDATASpecs(StructureSpecs):
    """Class defining the structure file format for LAMMPS data.

    Parameters
    ----------
    atom_style : {'full'}, optional
        LAMMPS atom style.

    """
    def __init__(self, atom_style='full'):
        super(LAMMPSDATASpecs, self).__init__()
        self._section_properties = OrderedDict()

        atoms_section_syntax = {}
        atoms_section_syntax['full'] = ['atomID', 'moleculeID', 'atomtype',
                                        'q', 'x', 'y', 'z', 'nx', 'ny', 'nz']
        atoms_section_syntax['atomic'] = ['atomID', 'moleculeID', 'atomtype',
                                          'x', 'y', 'z', 'nx', 'ny', 'nz']
        velocities_section_syntax = {}
        velocities_section_syntax['full'] = ['atomID', 'vx', 'vy', 'vz']

        syntax_dtypes = \
            {'atomID': int, 'moleculeID': int, 'atomtype': int, 'q': float,
             'mass': float, 'x': float, 'y': float, 'z': float,
             'nx': int, 'ny': int, 'nz': int,
             'vx': float, 'vy': float, 'vz': float}

        section_keys = ['Atoms', 'Masses', 'Velocities']
        self._section_syntax_dict = \
            {'Atoms': atoms_section_syntax[atom_style],
             'Masses': ['atomtype', 'mass'],
             'Velocities': velocities_section_syntax[atom_style]}

        for section_key in section_keys:
            self._section_properties[section_key] = OrderedDict()
            section_syntax_list = self._section_syntax_dict[section_key]
            for i, syntax in enumerate(section_syntax_list):
                self._section_properties[section_key][syntax] = \
                    {'dtype': syntax_dtypes[syntax],
                     'colnum': i+1,
                     'index': i}

        self._header_keys = ["atoms",
                             "atom types",
                             "bonds",
                             "bond types",
                             "angles",
                             "angle types",
                             "dihedrals",
                             "dihedral types",
                             "impropers",
                             "improper types",
                             "bodies",
                             "ellipsoids",
                             "lines",
                             "triangles",
                             "xlo xhi",
                             "ylo yhi",
                             "zlo zhi",
                             "xy xz yz"]

        self._headers = {'atoms': {'dtype': int, 'items': 1},
                         'atom types': {'dtype': int, 'items': 1},
                         'bonds': {'dtype': int, 'items': 1},
                         'bond types': {'dtype': int, 'items': 1},
                         'angles': {'dtype': int, 'items': 1},
                         'angle types': {'dtype': int, 'items': 1},
                         'dihedrals': {'dtype': int, 'items': 1},
                         'dihedral types': {'dtype': int, 'items': 1},
                         'impropers': {'dtype': int, 'items': 1},
                         'improper types': {'dtype': int, 'items': 1},
                         'bodies': {'dtype': int, 'items': 1},
                         'ellipsoids': {'dtype': int, 'items': 1},
                         'lines': {'dtype': int, 'items': 1},
                         'triangles': {'dtype': int, 'items': 1},
                         'xlo xhi': {'dtype': float, 'items': 2},
                         'ylo yhi': {'dtype': float, 'items': 2},
                         'zlo zhi': {'dtype': float, 'items': 2},
                         'xy xz yz': {'dtype': float, 'items': 3}}

        self._properties['headers'] = self._headers

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


class LAMMPSDUMPSpecs(StructureSpecs):
    """Class defining the structure file format for LAMMPS dump."""

    def __init__(self):
        pass
