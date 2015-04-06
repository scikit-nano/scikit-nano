# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.structures import UnrolledSWNT


def test1():
    unrolled_swnt = UnrolledSWNT(n=10, m=10)
    assert_equal(unrolled_swnt.n, 10)
    assert_equal(unrolled_swnt.m, 10)
    assert_equal(unrolled_swnt.element1, 'C')
    assert_equal(unrolled_swnt.element2, 'C')
    assert_equal(unrolled_swnt.nz, 1.0)
    assert_equal(unrolled_swnt.t1, 1)
    assert_equal(unrolled_swnt.t2, -1)
    assert_equal(unrolled_swnt.d, 10)
    assert_equal(unrolled_swnt.dR, 30)
    assert_equal(unrolled_swnt.N, 20)
    assert_equal(unrolled_swnt.R, (1, 0))
    assert_almost_equal(unrolled_swnt.chiral_angle, 30.0)
    assert_almost_equal(unrolled_swnt.Ch, 42.6, places=2)
    assert_almost_equal(unrolled_swnt.T, 2.46, places=2)


if __name__ == '__main__':
    nose.runmodule()
