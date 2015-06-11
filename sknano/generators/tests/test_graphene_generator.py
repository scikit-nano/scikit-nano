# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import GrapheneGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        generator = GrapheneGenerator(armchair_edge_length=5,
                                      zigzag_edge_length=5)
        generator.save()
        self.tmpdata.append(generator.fname)

    def test2(self):
        generator = GrapheneGenerator(armchair_edge_length=20,
                                      zigzag_edge_length=1)
        generator.save()
        self.tmpdata.append(generator.fname)

    def test3(self):
        generator = GrapheneGenerator(armchair_edge_length=1,
                                      zigzag_edge_length=20)
        generator.save()
        self.tmpdata.append(generator.fname)


if __name__ == '__main__':
    nose.runmodule()
