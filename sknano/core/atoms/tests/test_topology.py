#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal
import numpy as np
from sknano.generators import SWNTGenerator
from sknano.io import DATAReader
from sknano.core.math import Vector
from sknano.core.atoms import get_angle
from sknano.testing import AtomsTestFixture, generate_atoms


class TestCase(AtomsTestFixture):

    def test1(self):
        atoms = self.atoms
        atoms.update_attrs()

        atoms123 = atoms[:3]
        print(atoms123.data)
        angle = get_angle(atoms123.data)
        print(angle)

        # bonds = atoms.bonds
        # # angles = atoms.angles
        # assert_equal(bonds.Nbonds, atoms.coordination_numbers.sum())
        # print(bonds.mean_length)
        # print(bonds.atoms.Natoms)
        # print('bonds.mean_angle: {}'.format(bonds.mean_angle))


if __name__ == '__main__':
    nose.runmodule()
