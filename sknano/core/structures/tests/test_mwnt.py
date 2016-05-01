# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import warnings
warnings.simplefilter('always')

import nose
from nose.tools import assert_equal
from sknano.core.structures import MWNT


def test1():
    mwnt = MWNT(Ch=[(5, 5), (10, 10)])
    print(mwnt)
    print(mwnt.chiral_types)
    print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 2)
    assert_equal(mwnt.chiral_set, set(['armchair']))


def test2():
    mwnt = MWNT(max_walls=2)
    # print(mwnt)
    # print(mwnt.chiral_types)
    # print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 2)


def test3():
    mwnt = MWNT(Ch=[(5, 5), (10, 10), (15, 15), (20, 20)], Lz=10.0)
    # print(mwnt)
    # print(mwnt.chiral_types)
    # print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 4)
    assert_equal(mwnt.chiral_set, set(['armchair']))


def test4():
    mwnt = MWNT(Ch=[(5, 0), (5, 5), (10, 5), (10, 0), (10, 10)],
                Lz=5.0)
    # print(mwnt)
    # print(mwnt.chiral_types)
    # print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 5)
    assert_equal(mwnt.chiral_set, set(['armchair', 'chiral', 'zigzag']))


def test5():
    bundle = MWNT(Ch=[(5, 5), (10, 10)], nx=2, ny=2)
    assert_equal(bundle.Nwalls, 2)
    assert_equal(bundle.Ntubes, 4)


def test6():
    bundle = MWNT(max_walls=2, nx=5, ny=2)
    assert_equal(bundle.Nwalls, 2)
    assert_equal(bundle.Ntubes, 10)


def test7():
    Ch_list = [(5, 0), (10, 0)]
    mwnt = MWNT(Ch_list=Ch_list, Lz=10.0)
    print('\nmwnt.Natoms')
    print(mwnt.Natoms)
    # print('\nmwnt.Natoms_list')
    # print(mwnt.Natoms_list)
    print([swnt.unit_cell.basis.Natoms for swnt in mwnt.walls])
    # print('mwnt.walls')
    # print(mwnt.walls)
    print('\nmwnt.crystal_cell.Natoms')
    print(mwnt.crystal_cell.Natoms)

# def test8():
#     Ch_list = [(3, 3), (4, 4), (5, 5)]
#     mwnt = MWNT(Ch_list=Ch_list, Lz=10.0)


if __name__ == '__main__':
    nose.runmodule()
