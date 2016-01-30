# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal
# from sknano.io import XYZData, XYZReader, XYZWriter, XYZ2DATAConverter
from sknano.testing import IOTestFixture


class XYZTestFixture(IOTestFixture):

    def setUp(self):
        atoms = self.atoms = self.xyz_reader.atoms
        atoms.assign_unique_ids()
        atoms.assign_unique_types()
        atoms.update_attrs()


class Tests(XYZTestFixture):

    def test1(self):
        atoms = self.atoms
        assert_equal(atoms.Natoms, 40)
        assert_equal(list(set(atoms.atom_ids)), list(range(1, 41)))

    # def test2(self):
    #     # data = XYZData()
    #     # data.fpath = infile
    #     # data.read()
    #     # atoms = data.atoms
    #     atoms = self.atoms
    #     assert_not_equal(atoms.Natoms, 20)
    #     assert_equal(atoms.Natoms, 40)
    #     atoms.assign_unique_ids()
    #     assert_equal(list(set(atoms.atom_ids)), list(range(1, 41)))

    # def test3(self):
    #     atoms = self.atoms
    #     assert_not_equal(atoms.Natoms, 20)
    #     assert_equal(atoms.Natoms, 40)
    #     atoms.assign_unique_ids()
    #     assert_equal(list(set(atoms.atom_ids)), list(range(1, 41)))


if __name__ == '__main__':
    nose.runmodule()
