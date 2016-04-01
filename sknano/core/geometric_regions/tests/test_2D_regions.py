#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_false, assert_true, \
    assert_is_instance

import numpy as np

from sknano.core.math import Point, Points, Vector, Vectors
from sknano.testing import Geometric2DRegionsTestFixture


class Tests(Geometric2DRegionsTestFixture):

    def test_parallelogram(self):
        region0 = self.parallelogram
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.measure, 1)
        assert_is_instance(region0.centroid, Point)
        assert_equal(region0.centroid, Point([1, 0.5]))
        assert_true(region0.contains([0.5, 0.5]))
        assert_false(region0.contains([0, 1]))

        print('changing region0.o...')
        region0.o = [1, 1]
        # assert_equal(Vector(region0.u, p0=region0.o), region0.u)
        self.print_regions()

        print('rotating region0...')
        region0.rotate(angle=np.pi/2)
        assert_equal(region0.vectors, Vectors([region0.u, region0.v]))
        self.print_regions()

        print('changing region0.u...')
        region0.u = [1, 1]
        self.print_regions()

        assert_equal(region0.vectors, Vectors([region0.u, region0.v]))

        print('changing region0.o...')
        region0.o = [0, 0]
        self.print_regions()
        # assert_equal(Vector(region0.u, p0=region0.o), region0.u)

        print('rotating region0...')
        region0.rotate(angle=np.pi/2)
        self.print_regions()
        assert_equal(region0.vectors, Vectors([region0.u, region0.v]))

    def test_square(self):
        region0 = self.square
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.measure, 1)
        assert_equal(region0.centroid, Point([0, 0]))
        assert_true(region0.contains([0, 0]))
        assert_false(region0.contains([1, 1]))
        assert_is_instance(region0.centroid, Point)

    def test_rectangle(self):
        region0 = self.rectangle
        region_class = region0.__class__
        self.regions.append(region0)
        self.print_regions()

        assert_equal(region0.measure, 1)
        assert_equal(region0.centroid, Point([0.5, 0.5]))
        assert_true(region0.contains(region0.centroid))
        assert_false(region0.contains([1.1, 1.1]))
        assert_is_instance(region0.centroid, Point)

        pmin = [-10, -10]
        pmax = [10, 10]
        region1 = region_class(pmin=pmin, pmax=pmax)
        self.regions.append(region1)
        assert_true(np.allclose(region1.pmin, pmin))
        assert_true(np.allclose(region1.pmax, pmax))

        pmin = [5, 5]
        pmax = [10, 10]
        region2 = region_class(pmin=pmin, pmax=pmax)
        self.regions.append(region2)
        assert_true(np.allclose(region2.pmin, pmin))
        assert_true(np.allclose(region2.pmax, pmax))

        pmin = [-10, -10]
        pmax = [-5, -5]
        region3 = region_class(pmin=pmin, pmax=pmax)
        self.regions.append(region3)
        assert_true(np.allclose(region3.pmin, pmin))
        assert_true(np.allclose(region3.pmax, pmax))

    def test_circle(self):
        region0 = self.circle
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.measure, np.pi)
        assert_equal(region0.centroid, Point([0, 0]))
        assert_equal(region0.r, 1.0)
        assert_true(region0.contains([0, 0]))
        assert_false(region0.contains([1.1, 0]))
        assert_is_instance(region0.centroid, Point)

    def test_ellipse(self):
        region0 = self.ellipse
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.centroid, Point([0, 0]))
        assert_true(np.allclose(region0.measure, np.pi))
        assert_true(region0.contains(region0.centroid))
        assert_false(region0.contains([1.01, 1.01]))
        assert_is_instance(region0.centroid, Point)

    def test_triangle(self):
        region0 = self.triangle
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.centroid, Point([1 / 3, 1 / 3]))
        assert_true(np.allclose(region0.measure, 0.5))
        assert_equal(region0.p1, Point(nd=2))
        assert_equal(region0.p2, Point([0, 1]))
        assert_equal(region0.p3, Point([1, 0]))
        assert_true(region0.contains([0.25, 0.25]))
        assert_false(region0.contains([1, 1]))
        assert_is_instance(region0.centroid, Point)


if __name__ == '__main__':
    nose.runmodule()
