# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import BilayerGrapheneGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        generator = BilayerGrapheneGenerator(armchair_edge_length=5,
                                             zigzag_edge_length=5,
                                             layer_rotation_increment=45)
        generator.save()
        self.tmpdata.append(generator.fname)


if __name__ == '__main__':
    nose.runmodule()
