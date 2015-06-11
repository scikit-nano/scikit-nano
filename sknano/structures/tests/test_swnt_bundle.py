# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.structures import SWNTBundle


def test1():
    swnt = SWNTBundle(n=10, m=10)
    print(swnt)
    assert_equal(swnt.n, 10)
    assert_equal(swnt.m, 10)
    assert_equal(swnt.element1, 'C')
    assert_equal(swnt.element2, 'C')
    assert_equal(swnt.nz, 1.0)
    assert_equal(swnt.t1, 1)
    assert_equal(swnt.t2, -1)
    assert_equal(swnt.d, 10)
    assert_equal(swnt.dR, 30)
    assert_equal(swnt.N, 20)
    assert_equal(swnt.R, (1, 0))
    assert_almost_equal(swnt.chiral_angle, 30.0)
    assert_almost_equal(swnt.Ch, 42.6, places=2)
    assert_almost_equal(swnt.T, 2.46, places=2)
    assert_almost_equal(swnt.dt, 13.56, places=2)
    assert_almost_equal(swnt.rt, 6.78, places=2)
    assert_equal(swnt.Ntubes, 1)


def test2():
    bundle = SWNTBundle(n=10, m=10, nx=3, ny=3)
    print(bundle)
    assert_equal(bundle.element1, 'C')
    assert_equal(bundle.element2, 'C')
    assert_equal(bundle.Ntubes, 9)


if __name__ == '__main__':
    nose.runmodule()
