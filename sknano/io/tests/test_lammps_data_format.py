# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename
from sknano.io import DATAData, DATAReader, DATAWriter, DATA2XYZConverter, \
    DATAFormatSpec, atom_styles


def test_atom_styles():
    print('atom_styles: {}\n'.format(atom_styles))

    formatspec = DATAFormatSpec(atom_style='full')
    #print('formatspec.header_specs: {}\n'.format(header_specs))
    #print('formatspec.section_header_map: {}\n'.format(section_header_map))
    print('formatspec.section_attrs: {}\n'.format(formatspec.section_attrs))
    print('formatspec.section_attrs_specs: {}\n'.format(
        formatspec.section_attrs_specs))


def test_data():
    datafile = resource_filename('sknano', 'data/nanotubes/1010_1cell.data')
    atoms = DATAData(fpath=datafile).structure_data.atoms
    atoms.assign_unique_ids()
    atoms.update_attrs()
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)
    #assert_equal(DATAData(fpath=datafile).structure_data.atoms,
    #             DATAData(fpath=datafile).atoms)


def test_reader():
    datafile = resource_filename('sknano', 'data/nanotubes/1010_1cell.data')
    reader = DATAReader(datafile)
    atoms = reader.structure_data.atoms
    atoms.assign_unique_ids()
    atoms.update_attrs()
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)


if __name__ == '__main__':
    nose.runmodule()
