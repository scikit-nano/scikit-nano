# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import BilayerGrapheneGenerator
from sknano.testing import GeneratorTestFixture


class Tests(GeneratorTestFixture):

    def test1(self):
        generator = BilayerGrapheneGenerator(armchair_edge_length=1,
                                             zigzag_edge_length=1,
                                             layer_rotation_increment=45)
        generator.save()
        self.tmpdata.append(generator.fname)


if __name__ == '__main__':
    nose.runmodule()
