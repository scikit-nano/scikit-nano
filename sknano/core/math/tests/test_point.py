#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

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
    assert_equal(p.x, 1.0)
    assert_equal(p.y, 1.0)
    assert_equal(p.z, 1.0)

    p = Point([None, None, None])
    assert_true(np.allclose(p, np.zeros(3)))

    p = Point(np.arange(10))
    assert_true(np.allclose(p, np.arange(10)))


def test_rezero():
    p = Point([1e-9,1e-11,-1e-16])
    p.rezero()


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
    p.rotate(np.pi/2, rot_axis=[0, 0, 1], anchor_point=[1, 1, 0])
    assert_true(np.allclose(p, np.array([2, 1, 0])))


if __name__ == '__main__':
    nose.runmodule()
