#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

#import nose
import unittest

import numpy as np

from sknano.utils.geometric_shapes import Sphere, Ellipsoid


class TestSphere(unittest.TestCase):

    #def setUp(self):
    #    self.circle = Circle()

    def test_with_no_args(self):
        s = Sphere()
        print('sphere: {}'.format(s))
        print('sphere.center: {}'.format(s.center))
        self.assertTrue(np.allclose(s.center, np.zeros(3)))
        print('sphere.r: {}'.format(s.r))
        self.assertEqual(s.r, 1.0)

    def test_with_list_of_zeros(self):
        s = Sphere([0, 0, 0])
        print('sphere: {}'.format(s))
        print('sphere.center: {}'.format(s.center))
        print('sphere.r: {}'.format(s.r))
        self.assertTrue(np.allclose(s.center, np.zeros(3)))
        self.assertEqual(s.r, 1.0)

    def test_with_tuple_of_zeros(self):
        s = Sphere((0, 0, 0))
        print('sphere: {}'.format(s))
        print('sphere.center: {}'.format(s.center))
        print('sphere.r: {}'.format(s.r))
        self.assertTrue(np.allclose(s.center, np.zeros(3)))
        self.assertEqual(s.r, 1.0)


class TestEllipsoid(unittest.TestCase):

    def test_with_no_args(self):
        s = Ellipsoid()
        print('ellipsoid: {}'.format(s))
        print('ellipsoid.center: {}'.format(s.center))
        self.assertTrue(np.allclose(s.center, np.zeros(3)))

    def test_with_list_of_zeros(self):
        s = Ellipsoid([0, 0, 0])
        print('ellipsoid: {}'.format(s))
        print('ellipsoid.center: {}'.format(s.center))
        self.assertTrue(np.allclose(s.center, np.zeros(3)))



if __name__ == '__main__':
    unittest.main()
