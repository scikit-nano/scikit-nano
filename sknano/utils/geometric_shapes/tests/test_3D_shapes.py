#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
import numpy as np

from sknano.core.math import Vector
from sknano.utils.geometric_shapes import Parallelepiped, Sphere, Ellipsoid


def test_parallelepiped():
    s = Parallelepiped()
    assert_true(np.allclose(s.o, np.zeros(3)))

    assert_equal(s.volume, 1.0)


def test_sphere():
    s = Sphere()
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)

    s = Sphere([0, 0, 0])
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)

    s = Sphere((0, 0, 0))
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)


def test_ellipsoid():
    s = Ellipsoid()
    assert_true(np.allclose(s.center, np.zeros(3)))

    s = Ellipsoid([0, 0, 0])
    assert_true(np.allclose(s.center, np.zeros(3)))


if __name__ == '__main__':
    nose.runmodule()
