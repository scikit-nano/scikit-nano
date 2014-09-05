# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.structures import SWNT


def test1():
    swnt = SWNT(n=10, m=10)
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
    swnt = SWNT(n=20, m=10)
    assert_equal(swnt.element1, 'C')
    assert_equal(swnt.element2, 'C')
    assert_equal(swnt.n, 20)
    assert_equal(swnt.m, 10)
    assert_equal(swnt.nz, 1.0)
    assert_equal(swnt.t1, 4)
    assert_equal(swnt.t2, -5)
    assert_equal(swnt.d, 10)
    assert_equal(swnt.dR, 10)
    assert_equal(swnt.N, 140)
    assert_equal(swnt.R, (1, -1))
    assert_almost_equal(swnt.chiral_angle, 19.11, places=2)
    assert_almost_equal(swnt.Ch, 65.07, places=2)
    assert_almost_equal(swnt.T, 11.27, places=2)
    assert_almost_equal(swnt.dt, 20.71, places=2)
    assert_almost_equal(swnt.rt, 10.36, places=2)
    assert_equal(swnt.Ntubes, 1)


def test3():
    swnt = SWNT(n=20, m=0)
    assert_equal(swnt.element1, 'C')
    assert_equal(swnt.element2, 'C')
    assert_equal(swnt.n, 20)
    assert_equal(swnt.m, 0)
    assert_equal(swnt.nz, 1.0)
    assert_equal(swnt.t1, 1)
    assert_equal(swnt.t2, -2)
    assert_equal(swnt.d, 20)
    assert_equal(swnt.dR, 20)
    assert_equal(swnt.N, 40)
    assert_equal(swnt.R, (1, -1))
    assert_almost_equal(swnt.chiral_angle, 0.0, places=2)
    assert_almost_equal(swnt.Ch, 49.2, places=1)
    assert_almost_equal(swnt.T, 4.26, places=2)
    assert_almost_equal(swnt.dt, 15.7, places=1)
    assert_almost_equal(swnt.rt, 7.8, places=1)
    assert_equal(swnt.Ntubes, 1)


if __name__ == '__main__':
    nose.runmodule()
