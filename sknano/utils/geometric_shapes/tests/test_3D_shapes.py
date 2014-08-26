#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
import numpy as np

from sknano.utils.geometric_shapes import Sphere, Ellipsoid


def test_instantiation():
    s = Sphere()
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)

    s = Sphere([0, 0, 0])
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)

    s = Sphere((0, 0, 0))
    assert_true(np.allclose(s.center, np.zeros(3)))
    assert_equal(s.r, 1.0)

    s = Ellipsoid()
    assert_true(np.allclose(s.center, np.zeros(3)))

    s = Ellipsoid([0, 0, 0])
    assert_true(np.allclose(s.center, np.zeros(3)))


if __name__ == '__main__':
    nose.runmodule()
