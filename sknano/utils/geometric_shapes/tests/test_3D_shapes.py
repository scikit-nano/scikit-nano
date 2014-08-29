#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
import numpy as np

from sknano.utils.geometric_shapes import Parallelepiped, Cuboid, Cube, \
    Ellipsoid, Spheroid, Sphere


def test_parallelepiped():
    s = Parallelepiped()
    assert_is_instance(s, Parallelepiped)
    assert_true(np.allclose(s.o, np.zeros(3)))
    assert_equal(s.volume, 1.0)


def test_cuboid():
    s = Cuboid()
    assert_is_instance(s, Cuboid)


def test_cube():
    s = Cube()
    assert_is_instance(s, Cube)


def test_ellipsoid():
    s = Ellipsoid()
    assert_is_instance(s, Ellipsoid)
    assert_true(np.allclose(s.center, np.zeros(3)))

    s = Ellipsoid([0, 0, 0])
    assert_true(np.allclose(s.center, np.zeros(3)))


def test_spheroid():
    s = Spheroid()
    assert_is_instance(s, Spheroid)


def test_sphere():
    s = Sphere()
    assert_is_instance(s, Sphere)
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)

    s = Sphere([0, 0, 0])
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)

    s = Sphere((0, 0, 0))
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)


if __name__ == '__main__':
    nose.runmodule()
