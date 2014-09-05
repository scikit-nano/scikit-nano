# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.io import XYZData, XYZReader, XYZWriter, XYZ2DATAConverter


def test_reader():
    infile = '1010r_1cell.xyz'
    reader = XYZReader(infile)
    atoms = reader.atoms
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)
    assert_equal(set(atoms.atom_ids), {0})
    atoms.assign_unique_ids()
    assert_equal(list(set(atoms.atom_ids)), list(xrange(1, 41)))


def test_dataio1():
    infile = '1010r_1cell.xyz'
    data = XYZData()
    data.fpath = infile
    data.read()
    atoms = data.atoms
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)
    assert_equal(set(atoms.atom_ids), {0})
    atoms.assign_unique_ids()
    assert_equal(list(set(atoms.atom_ids)), list(xrange(1, 41)))


def test_dataio2():
    infile = '1010r_1cell.xyz'
    data = XYZData(infile)
    atoms = data.atoms
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)
    assert_equal(set(atoms.atom_ids), {0})
    atoms.assign_unique_ids()
    assert_equal(list(set(atoms.atom_ids)), list(xrange(1, 41)))


if __name__ == '__main__':
    nose.runmodule()
