# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_true, assert_is_instance
from sknano.generators import MWNTGenerator
from sknano.testing import GeneratorTestFixture


class Tests(GeneratorTestFixture):

    def test1(self):
        mwnt = MWNTGenerator(max_walls=2, Lz=1.0)
        assert_is_instance(mwnt.todict(), dict)
        assert_equal(mwnt.todict()['max_walls'], mwnt.Nwalls)
        mwnt.save()
        self.tmpdata.append(mwnt.fname)
        mwnt.save(structure_format='data')
        self.tmpdata.append(mwnt.fname)

    def test2(self):
        mwnt = MWNTGenerator(Ch=[(5, 5), (10, 10)], Lz=1.0)
        assert_equal(len(mwnt.todict()['Ch_list']), mwnt.Nwalls)
        assert_equal(mwnt.Nwalls, 2)
        assert_equal(mwnt.chiral_set, set(['armchair']))
        mwnt.save()
        self.tmpdata.append(mwnt.fname)
        mwnt.save(structure_format='data')
        self.tmpdata.append(mwnt.fname)

    def test3(self):
        mwnt = MWNTGenerator(Ch=[(3, 3), (5, 5), (10, 10)], Lz=1.0)
        assert_equal(mwnt.Nwalls, 3)
        assert_equal(len(mwnt.todict()['Ch_list']), mwnt.Nwalls)
        assert_equal(mwnt.chiral_set, set(['armchair']))
        mwnt.save()
        self.tmpdata.append(mwnt.fname)
        mwnt.save(structure_format='data')
        self.tmpdata.append(mwnt.fname)

    def test4(self):
        bundle = MWNTGenerator(max_walls=2, Lz=1.0, bundle_geometry='hexagon')
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test5(self):
        Ch = [(5, 5), (6, 6), (7, 7), (8, 8)]
        bundle = MWNTGenerator(Ch=Ch, Lz=0.5, bundle_geometry='hexagon')
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test6(self):
        Ch = [(3, 3), (4, 4), (5, 5)]
        bundle = MWNTGenerator(Ch=Ch, nx=5, ny=2, Lz=0.5)
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test7(self):
        Ch_list = [(3, 3), (4, 4), (5, 5)]
        mwnt = MWNTGenerator(Ch_list=Ch_list, Lz=1.0)
        bundle = MWNTGenerator(Ch_list=Ch_list, nx=3, ny=2, Lz=1.0)
        assert_true(bundle.Natoms, bundle.Ntubes * mwnt.Natoms)

    def test8(self):
        Ch_list = [(5, 0), (10, 0)]
        mwnt = MWNTGenerator(Ch_list=Ch_list, Lz=1.0, verbose=False)
        assert_true(mwnt.bounding_box.zmax < 11 and
                    mwnt.bounding_box.zmax > 10)
        mwnt.save()
        print('\nmwnt.Natoms')
        print(mwnt.Natoms)
        # print('\nmwnt.Natoms_list')
        # print(mwnt.Natoms_list)
        print([swnt.unit_cell.basis.Natoms for swnt in mwnt.walls])
        # print('mwnt.walls')
        # print(mwnt.walls)
        print('\nmwnt.crystal_cell.Natoms')
        print(mwnt.crystal_cell.Natoms)


if __name__ == '__main__':
    nose.runmodule()
