# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import MWNTGenerator


def test1():
    mwnt = MWNTGenerator(max_shells=3, Lz=1.0, fix_Lz=True)
    print(mwnt)
    print(mwnt.todict())
    mwnt.save_data()
    mwnt.save_data(structure_format='data')


def test2():
    mwnt = MWNTGenerator(Ch=[(10, 10), (50, 50)], Lz=1.0, fix_Lz=True)
    print(mwnt)
    print(mwnt.todict())
    assert_equal(mwnt.Nwalls, 2)
    assert_equal(mwnt.chiral_set, set(['armchair']))
    mwnt.save_data()
    mwnt.save_data(structure_format='data')


def test3():
    mwnt = MWNTGenerator(Ch=[(5, 5), (10, 10), (20, 20)], nz=5.0)
    print(mwnt)
    print(mwnt.todict())
    assert_equal(mwnt.Nwalls, 3)
    assert_equal(mwnt.chiral_set, set(['armchair']))
    mwnt.save_data()
    mwnt.save_data(structure_format='data')

if __name__ == '__main__':
    nose.runmodule()
