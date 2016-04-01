#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import unittest

import nose
from nose.tools import assert_equal, assert_false, assert_true, \
    assert_is_instance
import numpy as np

from sknano.core.math import Point, Vector, Vectors
from sknano.testing import Geometric3DRegionsTestFixture


class Tests(Geometric3DRegionsTestFixture):

    def test_parallelepiped1(self):
        region0 = self.parallelepiped
        region_class = region0.__class__
        region1 = \
            region_class(o=[1, 1, 0], u=[1, 0, 0], v=[1, 1, 0], w=[1, 1, 1])
        self.regions.extend([region0, region1])
        self.print_regions()
        assert_equal(region0.o, Point([0, 0, 0]))
        assert_true(np.allclose(region0.measure, 1.0))
        assert_equal(region0.centroid, Point([1, 1, 0.5]))
        assert_true(region0.contains([0.5, 0.5, 0.5]))
        assert_false(region0.contains([0, 0, 1]))
        assert_is_instance(region0.centroid, Point)

        print('changing region0.o...')
        region0.o = [1, 1, 1]
        self.print_regions()
        # assert_equal(Vector(region0.u, p0=region0.o), region0.u)

        print('rotating region0...')
        region0.rotate(angle=np.pi/2, axis='z')
        assert_equal(region0.vectors,
                     Vectors([region0.u, region0.v, region0.w]))
        self.print_regions()

        print('changing region0.u...')
        region0.u = [1, 1, 1]
        self.print_regions()
        assert_equal(region0.vectors,
                     Vectors([region0.u, region0.v, region0.w]))

        print('changing region0.o...')
        region0.o = [0, 0, 0]
        self.print_regions()
        # assert_equal(Vector(region0.u, p0=region0.o), region0.u)

        print('rotating r...')
        region0.rotate(angle=np.pi/2, axis='z')
        self.print_regions()

        assert_equal(region0.vectors,
                     Vectors([region0.u, region0.v, region0.w]))

    def test_cuboid(self):
        region0 = self.cuboid
        region_class = region0.__class__
        region1 = region_class(pmin=[-1, -1, -1], pmax=[1, 1, 1])
        self.regions.extend([region0, region1])
        self.print_regions()
        assert_equal(region0.centroid, Point([0.5, 0.5, 0.5]))
        assert_true(np.allclose(region0.measure, 1.0))
        assert_true(region0.contains(region0.centroid))
        assert_false(region0.contains([-0.1, 0, 0]))
        assert_is_instance(region0.centroid, Point)

        pmin = [-10, -10, -10]
        pmax = [10, 10, 10]
        region0 = region_class(pmin=pmin, pmax=pmax)
        assert_true(np.allclose(region0.pmin, pmin))
        assert_true(np.allclose(region0.pmax, pmax))

        pmin = [5, 5, 5]
        pmax = [10, 10, 10]
        region0 = region_class(pmin=pmin, pmax=pmax)
        assert_true(np.allclose(region0.pmin, pmin))
        assert_true(np.allclose(region0.pmax, pmax))

        pmin = [-10, -10, -10]
        pmax = [-5, -5, -5]
        region0 = region_class(pmin=pmin, pmax=pmax)
        assert_true(np.allclose(region0.pmin, pmin))
        assert_true(np.allclose(region0.pmax, pmax))

    def test_cube(self):
        region0 = self.cube
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.centroid, Point([0, 0, 0]))
        assert_true(np.allclose(region0.measure, 1.0))
        assert_true(region0.contains([0, 0, 0]))
        assert_false(region0.contains([1, 1, 1]))
        assert_is_instance(region0.centroid, Point)

    def test_ellipsoid(self):
        region0 = self.ellipsoid
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.centroid, Point([0, 0, 0]))
        assert_true(np.allclose(region0.measure, 4 / 3 * np.pi))
        assert_true(region0.contains(region0.centroid))
        assert_false(region0.contains([1.01, 1.01, 1.01]))
        assert_is_instance(region0.centroid, Point)

    def test_sphere(self):
        region0 = self.sphere
        region_class = region0.__class__
        region1 = region_class(center=[1, 1, 1], r=2.0)
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_equal(region0.center, Point([0, 0, 0]))
        assert_equal(region0.r, 1.0)
        assert_true(np.allclose(region0.measure, 4 / 3 * np.pi))
        assert_true(region0.contains([0, 0, 0]))
        assert_false(region0.contains([1.1, 1.1, 1.1]))
        assert_is_instance(region0.centroid, Point)

    def test_cylinder(self):
        region0 = self.cylinder
        region_class = region0.__class__
        region1 = region_class()

        self.regions.extend([region0, region1])
        self.print_regions()

        assert_true(np.allclose(region0.measure, 2 * np.pi))
        assert_equal(region0.p1, Point([0, 0, -1]))
        assert_equal(region0.p2, Point([0, 0, 1]))
        assert_true(region0.contains([0, 0, 0]))
        assert_false(region0.contains([1.1, 0, 0]))
        assert_equal(region0.centroid, Point())
        region0.p2 = [0, 0, 2]
        assert_equal(region0.centroid, Point([0, 0, 0.5]))
        assert_is_instance(region0.centroid, Point)
        assert_is_instance(region0.axis, Vector)

    def test_cone(self):
        region0 = self.cone
        region_class = region0.__class__
        region1 = region_class()
        self.regions.extend([region0, region1])
        self.print_regions()

        assert_true(np.allclose(region0.measure, 2 * np.pi / 3))
        assert_equal(region0.p1, Point(np.zeros(3)))
        assert_equal(region0.p2, Point([0, 0, 2]))
        assert_true(region0.contains([0, 0, 1]))
        assert_false(region0.contains([0, 0, 2.1]))
        assert_is_instance(region0.centroid, Point)
        assert_is_instance(region0.axis, Vector)


if __name__ == '__main__':
    nose.runmodule()
