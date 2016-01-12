# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal
from sknano.structures import MWNT


def test1():
    mwnt = MWNT(Ch=[(5, 5), (10, 10)])
    print(mwnt)
    print(mwnt.todict())
    print(mwnt.chiral_types)
    print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 2)
    assert_equal(mwnt.chiral_set, set(['armchair']))


def test2():
    mwnt = MWNT(max_walls=5)
    print(mwnt)
    print(mwnt.todict())
    print(mwnt.chiral_types)
    print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 5)


def test3():
    mwnt = MWNT(Ch=[(5, 5), (10, 10), (15, 15), (20, 20)], Lz=2.5)
    print(mwnt)
    print(mwnt.todict())
    print(mwnt.chiral_types)
    print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 4)
    assert_equal(mwnt.chiral_set, set(['armchair']))


def test4():
    mwnt = MWNT(Ch=[(5, 0), (5, 5), (10, 5), (10, 0), (10, 10), (20, 10)],
                Lz=1.0)
    print(mwnt)
    print(mwnt.todict())
    print(mwnt.chiral_types)
    print(mwnt.chiral_set)
    assert_equal(mwnt.Ntubes, 1)
    assert_equal(mwnt.Nwalls, 6)
    assert_equal(mwnt.chiral_set, set(['armchair', 'chiral', 'zigzag']))


def test5():
    bundle = MWNT(Ch=[(5, 5), (10, 10)], nx=2, ny=2)
    assert_equal(bundle.Nwalls, 2)
    assert_equal(bundle.Ntubes, 4)


def test6():
    bundle = MWNT(max_walls=5, nx=5, ny=2)
    assert_equal(bundle.Nwalls, 5)
    assert_equal(bundle.Ntubes, 10)


if __name__ == '__main__':
    nose.runmodule()
