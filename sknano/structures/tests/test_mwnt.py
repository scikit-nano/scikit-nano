# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.structures import MWNT


def test1():
    mwnt = MWNT(Ch=[(5,5),(10,10)])
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

if __name__ == '__main__':
    nose.runmodule()
