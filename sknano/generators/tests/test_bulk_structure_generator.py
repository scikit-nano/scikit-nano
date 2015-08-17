# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import assert_equal

import numpy as np

from sknano.core.crystallography import AlphaQuartz, DiamondStructure, \
    CaesiumChlorideStructure, RocksaltStructure, ZincblendeStructure, \
    Gold, Copper, MoS2
from sknano.generators import BCCStructureGenerator, FCCStructureGenerator, \
    AlphaQuartzGenerator, DiamondStructureGenerator, IronGenerator, \
    GoldGenerator, CopperGenerator, CaesiumChlorideStructureGenerator, \
    RocksaltStructureGenerator, ZincblendeStructureGenerator, MoS2Generator
from sknano.io import XYZReader, DATAReader
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        quartz = AlphaQuartzGenerator()
        quartz.save()
        self.tmpdata.append(quartz.fname)
        assert_equal(quartz.atoms.Natoms, AlphaQuartz().basis.Natoms)

    def test2(self):
        quartz = AlphaQuartzGenerator(scaling_matrix=[2, 2, 2])
        quartz.save()
        self.tmpdata.append(quartz.fname)
        assert_equal(quartz.atoms.Natoms, AlphaQuartz().basis.Natoms * 2 ** 3)

    def test3(self):
        gold = GoldGenerator()
        gold.save()
        self.tmpdata.append(gold.fname)
        assert_equal(gold.atoms.Natoms, Gold().basis.Natoms)

    def test4(self):
        gold = GoldGenerator(scaling_matrix=[2, 2, 2])
        gold.save()
        self.tmpdata.append(gold.fname)
        assert_equal(gold.atoms.Natoms, Gold().basis.Natoms * 2 ** 3)

    def test5(self):
        copper = CopperGenerator()
        copper.save()
        self.tmpdata.append(copper.fname)
        assert_equal(copper.atoms.Natoms, Copper().basis.Natoms)

    def test6(self):
        copper = CopperGenerator(scaling_matrix=[2, 2, 2])
        copper.save()
        self.tmpdata.append(copper.fname)
        assert_equal(copper.atoms.Natoms, Copper().basis.Natoms * 2 ** 3)

    def test7(self):
        diamond = DiamondStructureGenerator()
        diamond.save()
        self.tmpdata.append(diamond.fname)
        assert_equal(diamond.atoms.Natoms, DiamondStructure().basis.Natoms)

    def test8(self):
        diamond = DiamondStructureGenerator(scaling_matrix=[2, 2, 2])
        diamond.save()
        self.tmpdata.append(diamond.fname)
        assert_equal(diamond.atoms.Natoms,
                     DiamondStructure().basis.Natoms * 2 ** 3)

    def test9(self):
        caesium_chloride = CaesiumChlorideStructureGenerator()
        caesium_chloride.save()
        self.tmpdata.append(caesium_chloride.fname)
        assert_equal(caesium_chloride.atoms.Natoms,
                     CaesiumChlorideStructure().basis.Natoms)

    def test10(self):
        caesium_chloride = \
            CaesiumChlorideStructureGenerator(scaling_matrix=[2, 2, 2])
        caesium_chloride.save()
        self.tmpdata.append(caesium_chloride.fname)
        assert_equal(caesium_chloride.atoms.Natoms,
                     CaesiumChlorideStructure().basis.Natoms * 2 ** 3)

    def test11(self):
        molybdenum_disulphide = MoS2Generator()
        molybdenum_disulphide.save()
        self.tmpdata.append(molybdenum_disulphide.fname)
        assert_equal(molybdenum_disulphide.atoms.Natoms,
                     MoS2().basis.Natoms)

    def test12(self):
        molybdenum_disulphide = MoS2Generator(scaling_matrix=[2, 2, 2])
        molybdenum_disulphide.save()
        self.tmpdata.append(molybdenum_disulphide.fname)
        assert_equal(molybdenum_disulphide.atoms.Natoms,
                     MoS2().basis.Natoms * 2 ** 3)

    def test13(self):
        assert_equal(GoldGenerator().atoms, FCCStructureGenerator('Au').atoms)

    def test14(self):
        assert_equal(IronGenerator().atoms, BCCStructureGenerator('Fe').atoms)


if __name__ == '__main__':
    nose.runmodule()
