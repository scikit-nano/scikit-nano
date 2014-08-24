#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

#import nose
import unittest

import numpy as np

from sknano.core.math import Point, Vector, vector as vec


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


class TestVector(unittest.TestCase):
    def test_vector_with_no_args(self):
        v = Vector()
        self.assertTrue(np.allclose(v, np.zeros(3)))
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.zeros(3)))

    def test_vector_with_list_of_ones(self):
        v = Vector([1, 1, 1])
        self.assertTrue(np.allclose(v, np.ones(3)))
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 1.0)
        self.assertEqual(v.z, 1.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))

    def test_vector_with_kwargs(self):
        v = Vector(p=[1, 1, 1])
        self.assertTrue(np.allclose(v, np.ones(3)))
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 1.0)
        self.assertEqual(v.z, 1.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))

        v = Vector(p0=[1, 1, 1])
        self.assertTrue(np.allclose(v, -np.ones(3)))
        self.assertEqual(v.x, -1.0)
        self.assertEqual(v.y, -1.0)
        self.assertEqual(v.z, -1.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.zeros(3)))

        v = Vector(p=[1, 1, 1], p0=[1, 1, 1])
        self.assertTrue(np.allclose(v, np.zeros(3)))
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))

    def test_property_changes(self):
        v = Vector()
        v.p0 = np.ones(3)
        self.assertTrue(np.allclose(v, -np.ones(3)))
        self.assertEqual(v.x, -1.0)
        self.assertEqual(v.y, -1.0)
        self.assertEqual(v.z, -1.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.zeros(3)))

        v = Vector(p=[1., 1., 1.], p0=[1., 1., 1.])
        self.assertTrue(np.allclose(v, np.zeros(3)))
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))
        v.p0 = np.zeros(3)
        self.assertTrue(np.allclose(v, np.ones(3)))
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 1.0)
        self.assertEqual(v.z, 1.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))
        v.x = 5.0
        self.assertEqual(v.x, 5.0)
        self.assertEqual(v.x, v.p.x)
        v.y = -5.0
        self.assertEqual(v.y, -5.0)
        self.assertEqual(v.y, v.p.y)
        v.z = 0.0
        self.assertEqual(v.z, 0.0)
        self.assertEqual(v.z, v.p.z)

        v.p0 = np.array([0.5, -10.0, 2.5])
        self.assertEqual(v.p0.x, 0.5)
        self.assertEqual(v.p0.y, -10.0)
        self.assertEqual(v.p0.z, 2.5)

    def test_add(self):
        v1 = Vector([1.0, 0.0], p0=[5., 5.])
        v2 = Vector([1.0, 1.0], p0=[5., 5.])
        print('v1: {}'.format(v1))
        print('v2: {}\n'.format(v2))
        print('adding v1 + v2')
        v3 = v1 + v2
        print('v3 = v1 + v2: {}\n'.format(v3))
        self.assertIsInstance(v3, Vector)
        self.assertTrue(np.allclose(v3, np.array([2.0, 1.0])))
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        print('v3: {}\n'.format(v3))

        print('v3.nd: {}\n'.format(v3.nd))
        print('v3.p0: {}\n'.format(v3.p0))
        print('v3.p: {}\n'.format(v3.p))

        print('adding v2 + v1')
        v3 = v2 + v1
        print('v3 = v2 + v1: {}\n'.format(v3))
        self.assertIsInstance(v3, Vector)
        self.assertTrue(np.allclose(v3, np.array([2.0, 1.0])))
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        print('v3: {}\n'.format(v3))

        print('v3.nd: {}\n'.format(v3.nd))
        print('v3.p0: {}\n'.format(v3.p0))
        print('v3.p: {}\n'.format(v3.p))

    def test_2D_dot(self):
        v1 = Vector([1.0, 0.0])
        v2 = Vector([1.0, 1.0])

        # this way works, but defeats the purpose of the subclasses
        print('\ncomputing np.dot(np.asarray(v1), np.asarray(v2))')
        c = np.dot(np.asarray(v1), np.asarray(v2))
        print('c = np.dot(v1, v2): {}\n'.format(c))
        self.assertEqual(c, 1.0)

        # now try the sknano.core.math.vector.dot function
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        print('\ncomputing vec.dot(v1, v2)')
        c = vec.dot(v1, v2)
        print('c = vec.dot(v1, v2): {}\n'.format(c))
        self.assertEqual(c, 1.0)

    def test_2D_cross(self):
        v1 = Vector([1.0, 0.0])
        v2 = Vector([1.0, 1.0])

        # now try the sknano.core.math.vector.cross function
        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        self.assertEqual(v3, np.cross(np.asarray(v1), np.asarray(v2)))

        v1 = Vector([1, 0], p0=[1, 1])
        v2 = Vector([0, 1], p0=[1, 1])

        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        self.assertTrue(np.allclose(v3, np.cross(np.asarray(v1),
                                                 np.asarray(v2))))

        v1 = Vector([1, 0], p0=[1, 0])
        v2 = Vector([0, 1], p0=[0, 1])

        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        self.assertTrue(np.allclose(v3, np.cross(np.asarray(v1),
                                                 np.asarray(v2))))

    def test_3D_dot(self):
        v1 = Vector([1.0, 0.0, 0.0])
        v2 = Vector([1.0, 1.0, 0.0])

        # this way works, but defeats the purpose of the subclasses
        print('\ncomputing np.dot(np.asarray(v1), np.asarray(v2))')
        c = np.dot(np.asarray(v1), np.asarray(v2))
        print('c = np.dot(v1, v2): {}\n'.format(c))
        self.assertEqual(c, 1.0)

        # now try the sknano.core.math.vector.dot function
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        print('\ncomputing vec.dot(v1, v2)')
        c = vec.dot(v1, v2)
        print('c = vec.dot(v1, v2): {}\n'.format(c))
        self.assertEqual(c, 1.0)

    def test_3D_cross(self):
        v1 = Vector([1, 0, 0])
        v2 = Vector([0, 1, 0])

        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        print('v3.nd: {}\n'.format(v3.nd))
        print('v3.p0: {}\n'.format(v3.p0))
        print('v3.p: {}\n'.format(v3.p))
        self.assertIsInstance(v3, Vector)
        self.assertTrue(np.allclose(v3, np.cross(np.asarray(v1),
                                                 np.asarray(v2))))
        self.assertTrue(np.allclose(v3.p0, v1.p0))
        self.assertTrue(np.allclose(v3.p0, v2.p0))

        v1 = Vector([1, 0, 0], p0=[1, 1, 1])
        v2 = Vector([0, 1, 0], p0=[1, 1, 1])

        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        print('v3.nd: {}\n'.format(v3.nd))
        print('v3.p0: {}\n'.format(v3.p0))
        print('v3.p: {}\n'.format(v3.p))
        self.assertIsInstance(v3, Vector)
        self.assertTrue(np.allclose(v3, np.cross(np.asarray(v1),
                                                 np.asarray(v2))))
        self.assertTrue(np.allclose(v3.p0, v1.p0))
        self.assertTrue(np.allclose(v3.p0, v2.p0))

        v1 = Vector([1, 0, 0], p0=[1, 0, 0])
        v2 = Vector([0, 1, 0], p0=[0, 1, 0])

        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        print('v3.nd: {}\n'.format(v3.nd))
        print('v3.p0: {}\n'.format(v3.p0))
        print('v3.p: {}\n'.format(v3.p))
        self.assertIsInstance(v3, Vector)
        self.assertTrue(np.allclose(v3, np.cross(np.asarray(v1),
                                                 np.asarray(v2))))
        self.assertTrue(np.allclose(v3.p0, v1.p0))
        self.assertFalse(np.allclose(v3.p0, v2.p0))

        v1 = Vector([1, 2, 3], p0=[1, 0, 0])
        v2 = Vector([4, 5, 6], p0=[0, 1, 0])

        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        print('v3.nd: {}\n'.format(v3.nd))
        print('v3.p0: {}\n'.format(v3.p0))
        print('v3.p: {}\n'.format(v3.p))
        self.assertIsInstance(v3, Vector)
        self.assertTrue(np.allclose(v3, np.cross(np.asarray(v1),
                                                 np.asarray(v2))))
        self.assertTrue(np.allclose(v3.p0, v1.p0))
        self.assertFalse(np.allclose(v3.p0, v2.p0))

    def test_str(self):
        v1 = Vector()
        print('\n{!s}\n'.format(v1))

    def test_repr(self):
        v1 = Vector()
        print('\n{!r}\n'.format(v1))

    def test_transforms(self):
        v = Vector([1.0, 0.0, 0.0])
        v.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(v, np.array([0, 1, 0])))
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.array([0, 1, 0])))
        v.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(v, np.array([-1, 0, 0])))
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.array([-1, 0, 0])))
        v.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(v, np.array([0, -1, 0])))
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.array([0, -1, 0])))
        v.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(v, np.array([1, 0, 0])))
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.array([1, 0, 0])))

        v = Vector(p0=[1, 1, 1])
        self.assertTrue(np.allclose(v, -np.ones(3)))
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.array([0, 0, 0])))
        v.rotate(np.pi/2, rot_axis='z')
        self.assertTrue(np.allclose(v, np.array([1, -1, -1])))
        self.assertTrue(np.allclose(v.p0, np.array([-1, 1, 1])))
        self.assertTrue(np.allclose(v.p, np.zeros(3)))

        v = Vector(p0=[1, 1, 1])
        self.assertTrue(np.allclose(v, -np.ones(3)))
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.array([0, 0, 0])))
        print(v)
        v.rotate(np.pi/2, rot_axis='z', anchor_point=[1, 1, 1])
        print(v)
        self.assertTrue(np.allclose(v, np.array([1, -1, -1])))
        self.assertTrue(np.allclose(v.p0, np.array([1, 1, 1])))
        self.assertTrue(np.allclose(v.p, np.array([2, 0, 0])))

    def test_rotations(self):
        v = Vector(p0=[1, 1, 1])
        self.assertTrue(np.allclose(v, -np.ones(3)))
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.array([0, 0, 0])))
        v.rotate(np.pi/2, rot_axis='z', anchor_point=[1, 1, 1])
        self.assertTrue(np.allclose(v, np.array([1, -1, -1])))
        self.assertTrue(np.allclose(v.p0, np.array([1, 1, 1])))
        self.assertTrue(np.allclose(v.p, np.array([2, 0, 0])))

if __name__ == '__main__':
    unittest.main()
