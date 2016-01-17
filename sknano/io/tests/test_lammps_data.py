# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_not_equal
from pkg_resources import resource_filename
from sknano.io import DATAData, DATAReader, DATAFormatSpec, atom_styles


def test1():
    print('atom_styles: {}\n'.format(atom_styles))

    formatspec = DATAFormatSpec(atom_style='full')
    print('formatspec.section_attrs: {}\n'.format(formatspec.section_attrs))
    print('formatspec.section_attrs_specs: {}\n'.format(
        formatspec.section_attrs_specs))


def test2():
    datafile = resource_filename('sknano', 'data/nanotubes/1010_1cell.data')
    atoms = DATAData(fpath=datafile).atoms
    atoms.assign_unique_ids()
    atoms.update_attrs()
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)


def test3():
    datafile = resource_filename('sknano', 'data/nanotubes/1010_1cell.data')
    reader = DATAReader(datafile)
    atoms = reader.atoms
    atoms.assign_unique_ids()
    atoms.update_attrs()
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)


if __name__ == '__main__':
    nose.runmodule()
