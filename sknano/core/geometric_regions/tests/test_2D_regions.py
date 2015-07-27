#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

import numpy as np

from sknano.core.geometric_regions import Parallelogram, Rectangle, Square, \
    Ellipse, Circle, Triangle
from sknano.core.math import Point, Points, Vector, Vectors


def test_parallelogram():
    r = Parallelogram()
    assert_is_instance(r, Parallelogram)
    assert_equal(r.measure, 1)
    assert_equal(r.centroid, Point([1, 0.5]))
    assert_true(r.contains([0.5, 0.5]))
    assert_false(r.contains([0, 1]))
    r.rotate(angle=np.pi/2)
    assert_equal(r.vectors, Vectors([r.u, r.v]))
    r.u = [1, 1]
    assert_equal(r.vectors, Vectors([r.u, r.v]))

    print(r.u)
    r.o = [1, 1]
    print(r.u)
    assert_equal(Vector(r.u, p0=r.o), r.u)
    r.rotate(angle=np.pi/2)
    assert_equal(r.vectors, Vectors([r.u, r.v]))
    assert_is_instance(r.centroid, Point)


def test_square():
    r = Square()
    assert_is_instance(r, Square)
    assert_equal(r.measure, 1)
    assert_equal(r.centroid, Point([0, 0]))
    assert_true(r.contains([0, 0]))
    assert_false(r.contains([1, 1]))
    assert_is_instance(r.centroid, Point)


def test_rectangle():
    r = Rectangle()
    assert_is_instance(r, Rectangle)
    assert_equal(r.measure, 1)
    assert_equal(r.centroid, Point([0.5, 0.5]))
    assert_true(r.contains(r.centroid))
    assert_false(r.contains([1.1, 1.1]))
    assert_is_instance(r.centroid, Point)


def test_circle():
    r = Circle()
    assert_is_instance(r, Circle)
    assert_equal(r.measure, np.pi)
    assert_equal(r.centroid, Point([0, 0]))
    assert_equal(r.r, 1.0)
    assert_true(r.contains([0, 0]))
    assert_false(r.contains([1.1, 0]))
    assert_is_instance(r.centroid, Point)


def test_ellipse():
    r = Ellipse()
    assert_is_instance(r, Ellipse)
    assert_equal(r.centroid, Point([0, 0]))
    assert_true(np.allclose(r.measure, np.pi))
    assert_true(r.contains(r.centroid))
    assert_false(r.contains([1.01, 1.01]))
    assert_is_instance(r.centroid, Point)


def test_triangle():
    r = Triangle()
    assert_is_instance(r, Triangle)
    assert_equal(r.centroid, Point([1 / 3, 1 / 3]))
    assert_true(np.allclose(r.measure, 0.5))
    assert_equal(r.p1, Point(nd=2))
    assert_equal(r.p2, Point([0, 1]))
    assert_equal(r.p3, Point([1, 0]))
    assert_true(r.contains([0.25, 0.25]))
    assert_false(r.contains([1, 1]))
    assert_is_instance(r.centroid, Point)


if __name__ == '__main__':
    nose.runmodule()
