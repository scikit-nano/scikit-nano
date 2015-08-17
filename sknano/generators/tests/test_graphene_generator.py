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
        generator = GrapheneGenerator(armchair_edge_length=2,
                                      zigzag_edge_length=2)
        generator.save()
        self.tmpdata.append(generator.fname)

    def test2(self):
        generator = ConventionalCellGrapheneGenerator(armchair_edge_length=2,
                                                      zigzag_edge_length=2)
        generator.save(fname='2nmx2nm_1layer_conventional_cell_graphene-1.xyz')
        self.tmpdata.append(generator.fname)

    def test3(self):
        generator = \
            GrapheneGenerator.from_conventional_cell(armchair_edge_length=2,
                                                     zigzag_edge_length=2)
        generator.atoms.update_attrs()
        # generator.save(deepcopy=False)
        generator.save(fname='2nmx2nm_1layer_conventional_cell_graphene-2.xyz',
                       deepcopy=False)
        self.tmpdata.append(generator.fname)

    def test4(self):
        generator = PrimitiveCellGrapheneGenerator(edge_length=2)
        generator.save(fname='2nm_1layer_primitive_cell_graphene-1.xyz')
        self.tmpdata.append(generator.fname)
        # print(generator.Natoms)

    def test5(self):
        generator = \
            GrapheneGenerator.from_primitive_cell(edge_length=2)
        generator.save(fname='2nm_1layer_primitive_cell_graphene-2.xyz')
        self.tmpdata.append(generator.fname)

    def test6(self):
        generator = \
            GrapheneGenerator.from_conventional_cell(armchair_edge_length=2,
                                                     zigzag_edge_length=2,
                                                     nlayers=3)
        generator.save(fname='2nmx2nm_3layer_conventional_cell_graphene-3.xyz')
        self.tmpdata.append(generator.fname)

    def test7(self):
        generator = \
            GrapheneGenerator.from_primitive_cell(edge_length=2, nlayers=3)
        generator.save(fname='2nm_3layer_primitive_cell_graphene-3.xyz')
        self.tmpdata.append(generator.fname)

    def test8(self):
        generator = \
            GrapheneGenerator.from_conventional_cell(armchair_edge_length=2,
                                                     zigzag_edge_length=2,
                                                     nlayers=3)
        generator.save(fname='2nmx2nm_3layer_conventional_cell_graphene.data')
        self.tmpdata.append(generator.fname)

    def test9(self):
        generator = \
            GrapheneGenerator.from_primitive_cell(edge_length=2, nlayers=3)
        generator.save(fname='2nm_3layer_primitive_cell_graphene.data')
        self.tmpdata.append(generator.fname)


if __name__ == '__main__':
    nose.runmodule()
