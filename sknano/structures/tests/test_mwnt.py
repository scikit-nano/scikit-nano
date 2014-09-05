# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.structures import MWNT


def test1():
    mwnt = MWNT(Ch=[(5,5),(10,10)])
    assert_equal(mwnt.Nwalls, 2)

def test2():
    mwnt = MWNT(max_shells=5)
    assert_equal(mwnt.Nwalls, 5)
    assert_equal(mwnt.Ntubes, 1)

if __name__ == '__main__':
    nose.runmodule()
