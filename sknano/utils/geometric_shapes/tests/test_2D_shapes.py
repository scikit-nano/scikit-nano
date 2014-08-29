#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *

import numpy as np

from sknano.utils.geometric_shapes import Parallelogram, Rectangle, Square, \
    Ellipse, Circle


def test_parallelogram():
    s = Parallelogram()
    assert_is_instance(s, Parallelogram)


def test_square():
    s = Square()
    assert_is_instance(s, Square)


def test_rectangle():
    s = Rectangle()
    assert_is_instance(s, Rectangle)


def test_circle():
    c = Circle()
    assert_is_instance(c, Circle)
    assert_true(np.allclose(c.center, np.zeros(2)))
    assert_equal(c.r, 1.0)

    c = Circle([0, 0])
    assert_true(np.allclose(c.center, np.zeros(2)))
    assert_equal(c.r, 1.0)

    c = Circle((0, 0))
    assert_true(np.allclose(c.center, np.zeros(2)))
    assert_equal(c.r, 1.0)

    c = Circle((None, None))
    assert_true(np.allclose(c.center, np.zeros(2)))
    assert_equal(c.r, 1.0)


def test_ellipse():
    e = Ellipse()
    assert_is_instance(e, Ellipse)

    assert_true(np.allclose(e.center, np.zeros(2)))

    e = Ellipse([0, 0])
    assert_true(np.allclose(e.center, np.zeros(2)))

    e = Ellipse((0, 0))
    assert_true(np.allclose(e.center, np.zeros(2)))

    s = Rectangle()
    assert_true(np.allclose(s.center, np.zeros(2)))


if __name__ == '__main__':
    nose.runmodule()
