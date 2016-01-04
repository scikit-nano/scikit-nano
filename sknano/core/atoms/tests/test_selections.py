#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# from pkg_resources import resource_filename

import nose
from nose.tools import assert_equal, assert_equals, assert_false, \
    assert_true, assert_is_instance

import numpy as np

# import sknano.core.atoms
# from sknano.core.atoms import StructureAtom, StructureAtoms
# from sknano.generators import SWNTGenerator
# from sknano.io import DATAReader
# from sknano.structures import compute_Natoms
from sknano.core.atoms import SelectionParser
from sknano.testing import AtomsTestFixture


class TestCase(AtomsTestFixture):

    def test1(self):
        atoms = self.atoms
        atoms.kNN = 6
        atoms.NNrc = 9.0
        for atom in atoms:
            assert_equals(atom.CN, 0)
            atom.CN = 3
            assert_equals(atom.CN, 3)

        atoms.update_attrs()
        # selected_atoms = SelectionParser(atoms).parse("id 5 10 15")
        # print(selected_atoms)
        # selected_atoms = SelectionParser(atoms).parse("id 5 10 15 or z <=-20")
        # parser.parse("id 5 10 15 or type 1")
        # print(parser)
        selected_atoms = SelectionParser(atoms).parse("id 5 10 15 or z <=-20")
        print(selected_atoms.ids)
        # print(parser.ATOM_ATTRIBUTE)
        # print(atoms.select('(z >= -5) and (z <= 5)'))
        # atoms = atoms.filtered((atoms.z >= -5) & (atoms.z <= 5))
        # print('Natoms: {}'.format(atoms.Natoms))
        # for atom in atoms:
        #     assert_equals(atom.CN, atoms.kNN)

if __name__ == '__main__':
    nose.runmodule()
