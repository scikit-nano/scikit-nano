# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import GrapheneGenerator, \
    PrimitiveCellGrapheneGenerator, ConventionalCellGrapheneGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        generator = GrapheneGenerator(armchair_edge_length=5,
                                      zigzag_edge_length=5)
        generator.save()
        self.tmpdata.append(generator.fname)

    def test2(self):
        generator = ConventionalCellGrapheneGenerator(armchair_edge_length=5,
                                                      zigzag_edge_length=5)
        generator.save()
        self.tmpdata.append(generator.fname)

    def test3(self):
        generator = \
            GrapheneGenerator.from_conventional_cell(armchair_edge_length=5,
                                                     zigzag_edge_length=5)
        generator.atoms.update_attrs()
        generator.save(deepcopy=False)
        self.tmpdata.append(generator.fname)

    def test4(self):
        generator = PrimitiveCellGrapheneGenerator(edge_length=5)
        generator.save()
        self.tmpdata.append(generator.fname)
        print(generator.Natoms)

    def test5(self):
        generator = \
            GrapheneGenerator.from_primitive_cell(edge_length=5)
        generator.save()
        self.tmpdata.append(generator.fname)


if __name__ == '__main__':
    nose.runmodule()
