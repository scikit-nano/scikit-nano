#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
import numpy as np
import unittest

from sknano.core.math import rotation_matrix, transformation_matrix, \
    rotation_transform, Rx, Ry, Rz, Point, Vector


class TestTransforms(unittest.TestCase):

    def test_rotation_matrix(self):
        self.assertTrue(np.allclose(rotation_matrix(np.pi/2),
                                    np.array([[0.0, -1.0], [1.0, 0.0]])))

        sqrt2 = np.sqrt(2)
        self.assertTrue(np.allclose(rotation_matrix(np.pi/4),
                                    np.array([[1 / sqrt2, -1 / sqrt2],
                                              [1 / sqrt2, 1 / sqrt2]])))

        sqrt3 = np.sqrt(3)
        self.assertTrue(np.allclose(
            rotation_matrix(np.pi/2, rot_axis=Vector(np.ones(3))),
            np.array([[1 / 3, 1 / 3 - 1 / sqrt3, 1 / 3 + 1 / sqrt3],
                      [1 / 3 + 1 / sqrt3, 1 / 3, 1 / 3 - 1 / sqrt3],
                      [1 / 3 - 1 / sqrt3, 1 / 3 + 1 / sqrt3, 1 / 3]])))

        [[self.assertTrue(np.allclose(
            Rmat(angle), rotation_matrix(angle, rot_axis=axis)))
            for angle in np.linspace(0.0, 2 * np.pi, np.pi / 4)]
            for Rmat, axis in zip((Rx, Ry, Rz), ('x', 'y', 'z'))]

    def test_transformation_matrix(self):
        self.assertTrue(np.allclose(
            transformation_matrix(np.pi/2, rot_axis=Point([1.0, 1.0])),
            np.array([[0.0, -1.0, 2.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])))

    def test_point_rotation(self):
        self.assertTrue(np.allclose(
            rotation_transform(Point([1.0, 0.0]), np.pi / 2),
            np.array([0.0, 1.0])))

        with self.assertRaises(TypeError):
            np.allclose(rotation_transform(
                Point([1.0, 0.0, 0.0]), angle=np.pi/2),
                np.array([0.0, 1.0, 0.0]))

        self.assertTrue(np.allclose(
            rotation_transform(Point([1.0, 0.0, 0.0]), angle=np.pi / 2,
                               rot_axis=Vector([1.0, 0.0, 0.0])),
            np.array([1.0, 0.0, 0.0])))

        self.assertTrue(np.allclose(
            rotation_transform(Point([0.0, 1.0, 0.0]), angle=np.pi / 2,
                               rot_axis=Vector([1.0, 0.0, 0.0])),
            np.array([0.0, 0.0, 1.0])))


if __name__ == '__main__':
    unittest.main()
