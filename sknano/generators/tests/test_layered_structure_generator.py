#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import os

import nose
from nose.tools import assert_equal

# import numpy as np

from sknano.testing import GeneratorTestFixture
from sknano.generators import LayeredStructureGenerator


class Tests(GeneratorTestFixture):

    def test1(self):
        config = os.path.join(os.path.dirname(__file__),
                              'test_layered_structure.cfg')
        generator = LayeredStructureGenerator(config)
        generator.save()
        self.tmpdata.append(generator.fname)
        structure = generator.structure
        silver = structure.select('element Ag')
        carbon = structure.select('element C')
        assert_equal(silver.Natoms + carbon.Natoms, structure.Natoms)
        assert_equal(structure.Natoms, structure.select('all').Natoms)

if __name__ == '__main__':
    nose.runmodule()
