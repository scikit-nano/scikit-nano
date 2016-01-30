# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from collections import OrderedDict

import nose
from nose.tools import assert_equal, assert_not_equal, assert_is_instance
from sknano.testing import IOTestFixture

from sknano.io import DATAFormatSpec, atom_styles


class DATATestFixture(IOTestFixture):

    def setUp(self):
        atoms = self.atoms = self.data_reader.atoms
        atoms.assign_unique_ids()
        atoms.assign_unique_types()
        atoms.update_attrs()


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


if __name__ == '__main__':
    nose.runmodule()
