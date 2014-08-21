#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

#import nose
import unittest

import numpy as np

from sknano.core.npmathobj import Point, Vector, vector as vec


class TestPoint(unittest.TestCase):

    #def setUp(self):
    #    self.nppoint = NPPoint()

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


class TestVector(unittest.TestCase):
    def test_vector_with_no_args(self):
        v = Vector(verbose=True)
        print('vector: {}\n'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}\n'.format(v.p))
        self.assertTrue(np.allclose(v, np.zeros(3)))
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.zeros(3)))

    def test_vector_with_list_of_ones(self):
        v = Vector([1, 1, 1], verbose=True)
        print('vector: {}\n'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}\n'.format(v.p))
        self.assertTrue(np.allclose(v, np.ones(3)))
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 1.0)
        self.assertEqual(v.z, 1.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, v))

    def test_vector_with_p_set(self):
        v = Vector(p=[1, 1, 1], verbose=True)
        print('vector: {}\n'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}\n'.format(v.p))
        self.assertTrue(np.allclose(v, np.ones(3)))
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 1.0)
        self.assertEqual(v.z, 1.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, v))

    def test_vector_with_p0_set(self):
        v = Vector(p0=[1, 1, 1], verbose=True)
        print('vector: {}\n'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}\n'.format(v.p))
        self.assertTrue(np.allclose(v, -np.ones(3)))
        self.assertEqual(v.x, -1.0)
        self.assertEqual(v.y, -1.0)
        self.assertEqual(v.z, -1.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.zeros(3)))

    def test_vector_with_p_and_p0_set(self):
        v = Vector(p=[1, 1, 1], p0=[1, 1, 1])
        print('vector: {}\n'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}\n'.format(v.p))
        self.assertTrue(np.allclose(v, np.zeros(3)))
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))

    def test_changing_p0(self):
        v = Vector(verbose=True)
        print('vector: {}\n'.format(v))
        v.p0 = np.ones(3)
        print('vector: {}'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}\n'.format(v.p))
        self.assertTrue(np.allclose(v, -np.ones(3)))
        self.assertEqual(v.x, -1.0)
        self.assertEqual(v.y, -1.0)
        self.assertEqual(v.z, -1.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.zeros(3)))

    def test_property_changes(self):
        v = Vector(p=[1., 1., 1.], p0=[1., 1., 1.], verbose=True)
        print('vector: {}\n'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}\n'.format(v.p))
        self.assertTrue(np.allclose(v, np.zeros(3)))
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)
        self.assertTrue(np.allclose(v.p0, np.ones(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))
        print('changing v.p0')
        v.p0 = np.zeros(3)
        print('vector: {}'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}'.format(v.p))
        self.assertTrue(np.allclose(v, np.ones(3)))
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 1.0)
        self.assertEqual(v.z, 1.0)
        self.assertTrue(np.allclose(v.p0, np.zeros(3)))
        self.assertTrue(np.allclose(v.p, np.ones(3)))
        print('changing v.x')
        v.x = 5.0
        print('changing v.y')
        v.y = -5.0
        print('changing v.z')
        v.z = 0.0
        self.assertEqual(v.x, 5.0)
        self.assertEqual(v.y, -5.0)
        self.assertEqual(v.z, 0.0)
        print('vector: {}'.format(v))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print('vector.p0: {}'.format(v.p0))
        print('vector.p: {}'.format(v.p))
        print('vector.p.x: {}'.format(v.p.x))
        print('vector.p.y: {}'.format(v.p.y))
        print('vector.p.z: {}\n'.format(v.p.z))
        self.assertEqual(v.x, v.p.x)
        self.assertEqual(v.y, v.p.y)
        self.assertEqual(v.z, v.p.z)
        print('changing v.p0')
        v.p0 = np.array([0.5,-10,2.5])
        print('vector.p: {}'.format(v.p))
        print('vector.p0: {}'.format(v.p0))
        print('vector v=p-p0: {}'.format(v))
        print('vector.p.x: {}'.format(v.p.x))
        print('vector.p.y: {}'.format(v.p.y))
        print('vector.p.z: {}'.format(v.p.z))
        print('vector.p0.x: {}'.format(v.p0.x))
        print('vector.p0.y: {}'.format(v.p0.y))
        print('vector.p0.z: {}'.format(v.p0.z))
        print('vector.x: {}'.format(v.x))
        print('vector.y: {}'.format(v.y))
        print('vector.z: {}'.format(v.z))
        print()

    def test_add(self):
        v1 = Vector([1.0, 0.0], p0=[5., 5.], verbose=True)
        v2 = Vector([1.0, 1.0], p0=[5., 5.], verbose=True)
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

    def test_cross(self):
        v1 = Vector([1.0, 0.0], verbose=True)
        v2 = Vector([1.0, 1.0], verbose=True)

        # now try the sknano.core.npmathobj.vector.cross function
        print('\ncomputing vec.cross(v1, v2)')
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        v3 = vec.cross(v1, v2)
        print('v3 = vec.cross(v1, v2) = {}\n'.format(v3))
        self.assertIsInstance(v3, np.number)
        self.assertEqual(v3, np.cross(np.asarray(v1), np.asarray(v2)))

        v1 = Vector([1, 0, 0], p0=[1, 1, 1], verbose=True)
        v2 = Vector([0, 1, 0], p0=[1, 1, 1], verbose=True)

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

        v1 = Vector([1, 0, 0], p0=[1, 0, 0], verbose=True)
        v2 = Vector([0, 1, 0], p0=[0, 1, 0], verbose=True)

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

        v1 = Vector([1, 2, 3], p0=[1, 0, 0], verbose=True)
        v2 = Vector([4, 5, 6], p0=[0, 1, 0], verbose=True)

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

    def test_dot(self):
        v1 = Vector([1.0, 0.0], verbose=True)
        v2 = Vector([1.0, 1.0], verbose=True)

        # this way works, but defeats the purpose of the subclasses
        print('\ncomputing np.dot(np.asarray(v1), np.asarray(v2))')
        c = np.dot(np.asarray(v1), np.asarray(v2))
        print('c = np.dot(v1, v2): {}\n'.format(c))
        self.assertEqual(c, 1.0)

        # now try the sknano.core.npmathobj.vector.dot function
        print('v1: {}'.format(v1))
        print('v2: {}'.format(v2))
        print('\ncomputing vec.dot(v1, v2)')
        c = vec.dot(v1, v2)
        print('c = vec.dot(v1, v2): {}\n'.format(c))
        self.assertEqual(c, 1.0)

    def test_str(self):
        v1 = Vector()
        print('\n{!s}\n'.format(v1))

    def test_repr(self):
        v1 = Vector()
        print('\n{!r}\n'.format(v1))

if __name__ == '__main__':
    unittest.main()
