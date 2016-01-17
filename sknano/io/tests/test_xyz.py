# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_not_equal
from pkg_resources import resource_filename
from sknano.io import XYZData, XYZReader  # , XYZWriter, XYZ2DATAConverter


def test1():
    infile = resource_filename('sknano', 'data/nanotubes/1010_1cell.xyz')
    reader = XYZReader(infile)
    atoms = reader.atoms
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)
    atoms.assign_unique_ids()
    assert_equal(list(set(atoms.atom_ids)), list(range(1, 41)))


def test2():
    infile = resource_filename('sknano', 'data/nanotubes/1010_1cell.xyz')
    data = XYZData()
    data.fpath = infile
    data.read()
    atoms = data.atoms
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)
    atoms.assign_unique_ids()
    assert_equal(list(set(atoms.atom_ids)), list(range(1, 41)))


def test3():
    infile = resource_filename('sknano', 'data/nanotubes/1010_1cell.xyz')
    data = XYZData(infile)
    atoms = data.atoms
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)
    atoms.assign_unique_ids()
    assert_equal(list(set(atoms.atom_ids)), list(range(1, 41)))


if __name__ == '__main__':
    nose.runmodule()
