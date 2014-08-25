#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

#import nose
import unittest

import numpy as np

from sknano.utils.geometric_shapes import Parallelogram, Rectangle, Square, \
    Ellipse, Circle


class TestParallelogram(unittest.TestCase):
    def test_no_args(self):
        s = Parallelogram()
        print('parallelogram: {}'.format(s))


class TestCircle(unittest.TestCase):

    #def setUp(self):
    #    self.circle = Circle()

    def test_with_no_args(self):
        c = Circle()
        print('circle: {}'.format(c))
        print('circle.center: {}'.format(c.center))
        self.assertTrue(np.allclose(c.center, np.zeros(2)))
        print('circle.r: {}'.format(c.r))
        self.assertEqual(c.r, 1.0)

    def test_with_list_of_zeros(self):
        c = Circle([0, 0])
        print('circle: {}'.format(c))
        print('circle.center: {}'.format(c.center))
        self.assertTrue(np.allclose(c.center, np.zeros(2)))
        print('circle.r: {}'.format(c.r))
        self.assertEqual(c.r, 1.0)

    def test_with_tuple_of_zeros(self):
        c = Circle((0, 0))
        print('circle: {}'.format(c))
        print('circle.center: {}'.format(c.center))
        self.assertTrue(np.allclose(c.center, np.zeros(2)))
        print('circle.r: {}'.format(c.r))
        self.assertEqual(c.r, 1.0)

    def test_with_tuple_of_Nones(self):
        c = Circle((None, None))
        print('circle: {}'.format(c))
        print('circle.center: {}'.format(c.center))
        self.assertTrue(np.allclose(c.center, np.zeros(2)))
        print('circle.r: {}'.format(c.r))
        self.assertEqual(c.r, 1.0)


class TestEllipse(unittest.TestCase):

    def test_with_no_args(self):
        e = Ellipse()
        print('ellipse: {}'.format(e))
        print('ellipse.center: {}'.format(e.center))
        self.assertTrue(np.allclose(e.center, np.zeros(2)))

    def test_with_list_of_zeros(self):
        e = Ellipse([0, 0])
        print('ellipse: {}'.format(e))
        print('ellipse.center: {}'.format(e.center))
        self.assertTrue(np.allclose(e.center, np.zeros(2)))

    def test_with_tuple_of_zeros(self):
        e = Ellipse((0, 0))
        print('ellipse: {}'.format(e))
        print('ellipse.center: {}'.format(e.center))
        self.assertTrue(np.allclose(e.center, np.zeros(2)))


class TestRectangle(unittest.TestCase):

    def test_with_no_args(self):
        s = Rectangle()
        print('rectangle: {}'.format(s))
        print('rectangle.center: {}'.format(s.center))
        self.assertTrue(np.allclose(s.center, np.zeros(2)))


if __name__ == '__main__':
    unittest.main()
