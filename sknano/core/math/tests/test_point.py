#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

import numpy as np

from sknano.core.math import Point


def test_instantiation():
    p = Point()
    assert_true(np.allclose(p, np.zeros(3)))

    p = Point([0,0,0], dtype=float)
    assert_true(np.allclose(p, np.zeros(3)))

    p = Point(np.ones(3))
    assert_true(np.allclose(p, np.ones(3)))

    p = Point([None, None, None])
    assert_true(np.allclose(p, np.zeros(3)))

    p = Point(np.arange(10))
    assert_true(np.allclose(p, np.arange(10)))


def test_rezero():
    p = Point([1e-9,1e-11,-1e-16])
    p.rezero(epsilon=1e-10)
    assert_not_equal(p.x, 0.0)
    assert_equal(p.y, 0.0)
    assert_equal(p.z, 0.0)


def test_transforms():
    p = Point([1.0, 0.0, 0.0])
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([0, 1, 0])))
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([-1, 0, 0])))
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([0, -1, 0])))
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([1, 0, 0])))
    p.rotate(np.pi/2, rot_axis=[0, 0, 1], rot_point=[1, 1, 0])
    assert_true(np.allclose(p, np.array([2, 1, 0])))


def test_total_ordering():
    p0 = Point(np.zeros(3))
    p1 = Point(np.ones(3))
    p2 = Point(2 * np.ones(3))
    pm2 = Point(-2 * np.ones(3))
    plist = [p2, p1, p0, pm2]
    assert_equal(plist.index(p0), 2)
    assert_true(p0 == Point(p0))
    assert_equal(plist.index(p1), 1)
    assert_true(p1 == Point(p1))
    assert_equal(plist.index(p2), 0)
    assert_true(p2 == Point(p2))
    assert_equal(plist.index(pm2), 3)
    assert_true(pm2 == Point(pm2))
    plist.sort()
    assert_equal(plist.index(p0), 0)
    assert_true(p0 == plist[0])
    assert_equal(plist.index(p1), 1)
    assert_true(p1 == plist[1])
    assert_equal(plist.index(p2), 2)
    assert_true(p2 == plist[2])
    assert_equal(plist.index(pm2), 3)
    assert_true(pm2 == plist[-1])


if __name__ == '__main__':
    nose.runmodule()
