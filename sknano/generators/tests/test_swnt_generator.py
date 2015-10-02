# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import *

import numpy as np

from sknano.generators import SWNTGenerator
from sknano.io import XYZReader, DATAReader
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        swnt = SWNTGenerator(n=10, m=10)
        swnt.save()
        self.tmpdata.append(swnt.fname)
        swnt.save(structure_format='data')
        self.tmpdata.append(swnt.fname)

    def test2(self):
        swnt = SWNTGenerator(n=10, m=10, Lz=1.0, fix_Lz=True)
        swnt.save()
        self.tmpdata.append(swnt.fname)
        swnt.save(structure_format='data')
        self.tmpdata.append(swnt.fname)

    def test3(self):

        new_swnt = SWNTGenerator((10, 10))
        new_swnt.atoms.center_centroid()
        test_swnt = \
            XYZReader(resource_filename('sknano',
                                        'data/nanotubes/1010_1cell.xyz'))
        print(new_swnt.atoms.coords)
        print(test_swnt.atoms.coords)
        assert_true(np.allclose(new_swnt.atoms.coords, test_swnt.atoms.coords))


if __name__ == '__main__':
    nose.runmodule()
