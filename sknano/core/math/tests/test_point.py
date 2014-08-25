#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

#import nose
import unittest

import numpy as np

from sknano.core.math import Point


class TestPoint(unittest.TestCase):

    def test_point_with_no_args(self):
        pt = Point()
        print('point: {}'.format(pt))
        print('point.x: {}'.format(pt.x))
        print('point.y: {}'.format(pt.y))
        print('point.z: {}\n'.format(pt.z))

    def test_point_with_list_of_zeros(self):
        pt = Point([0,0,0], dtype=float)
        print('point: {}'.format(pt))
        print('point.x: {}'.format(pt.x))
        print('point.y: {}'.format(pt.y))
        print('point.z: {}\n'.format(pt.z))
        self.assertTrue(np.allclose(pt, np.zeros(3)))
        self.assertEqual(pt.x, 0.0)
        self.assertEqual(pt.y, 0.0)
        self.assertEqual(pt.z, 0.0)

    def test_point_with_ones(self):
        pt = Point(np.ones(3))
        print('point: {}'.format(pt))
        print('point.x: {}'.format(pt.x))
        print('point.y: {}'.format(pt.y))
        print('point.z: {}\n'.format(pt.z))
        self.assertTrue(np.allclose(pt, np.ones(3)))
        self.assertEqual(pt.x, 1.0)
        self.assertEqual(pt.y, 1.0)
        self.assertEqual(pt.z, 1.0)

    def test_point_with_list_of_None(self):
        pt = Point([None, None, None])
        print('point: {}'.format(pt))
        pt.x = 2.5
        print('point: {}'.format(pt))

    def test_rezero_coords(self):
        pt = Point([1e-9,1e-11,-1e-16])
        print('point: {}'.format(pt))
        print('point.x: {}'.format(pt.x))
        print('point.y: {}'.format(pt.y))
        print('point.z: {}'.format(pt.z))
        pt.rezero_coords()
        print('point: {}\n'.format(pt))

    def test_point_with_arange10(self):
        pt = Point(np.arange(10))
        print('point: {}\n'.format(pt))

    def test_transforms(self):
        p = Point([1.0, 0.0, 0.0])
        p.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(p, np.array([0, 1, 0])))
        p.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(p, np.array([-1, 0, 0])))
        p.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(p, np.array([0, -1, 0])))
        p.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(p, np.array([1, 0, 0])))
        p.rotate(np.pi/2, rot_axis=[0, 0, 1], anchor_point=[1, 1, 0])
        self.assertTrue(np.allclose(p, np.array([2, 1, 0])))


if __name__ == '__main__':
    unittest.main()
