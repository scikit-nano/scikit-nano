# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from collections import OrderedDict

import nose
from nose.tools import assert_equal, assert_not_equal, assert_is_instance
from sknano.testing import IOTestFixture, generate_structure

from sknano.io import DATAData, DATAReader, DATAWriter, DATAFormatSpec, \
    atom_styles


class DATATestFixture(IOTestFixture):

    @property
    def atoms(self):
        atoms = self.data_reader.atoms
        atoms.assign_unique_ids()
        atoms.assign_unique_types()
        atoms.update_attrs()
        return atoms


class Tests(DATATestFixture):

    def test1(self):
        print('atom_styles: {}\n'.format(atom_styles))
        formatspec = DATAFormatSpec(atom_style='full')
        assert_is_instance(formatspec.section_attrs, OrderedDict)
        assert_is_instance(formatspec.section_attrs_specs, OrderedDict)

    def test2(self):
        atoms = self.atoms
        assert_not_equal(atoms.Natoms, 20)
        assert_equal(atoms.Natoms, 40)

    def test3(self):
        data = self.datadata
        # print(data)
        # print(data.headers)
        # print(data.sections)
        atoms = data.atoms
        testfile1 = 'test3-1.data'
        self.tmpdata.append(testfile1)
        data.write(testfile1)
        test_atoms1 = DATAData(testfile1).atoms
        assert_equal(atoms, test_atoms1)
        testfile2 = 'test3-2.data'
        self.tmpdata.append(testfile2)
        DATAWriter.write(datafile=testfile2, atoms=atoms)
        test_atoms2 = DATAReader(testfile2).atoms
        assert_equal(atoms, test_atoms2)
        testfile3 = 'test3-3.data'
        self.tmpdata.append(testfile3)
        swnt = generate_structure(generator_class='SWNTGenerator',
                                  Ch=(10, 10))
        DATAWriter.write(testfile3, structure=swnt.structure,
                         atoms=swnt.atoms)
        test_atoms3 = DATAReader(testfile3).atoms
        # print(test_atoms3.ids)
        assert_equal(atoms, test_atoms3)


if __name__ == '__main__':
    nose.runmodule()
