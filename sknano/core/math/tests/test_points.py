#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_true, assert_equal, assert_not_equal, \
    assert_is_instance

import numpy as np

from sknano.core.math import Point, Points


def test1():
    p = Point()
    assert_true(np.allclose(p, np.zeros(3)))

    p = Point([0, 0, 0], dtype=float)
    assert_true(np.allclose(p, np.zeros(3)))

    p = Point(np.ones(3))
    assert_true(np.allclose(p, np.ones(3)))

    p = Point([None, None, None])
    assert_true(np.allclose(p, np.zeros(3)))

    p = Point(np.arange(10))
    assert_true(np.allclose(p, np.arange(10)))


def test2():
    p = Point([1e-9, 1e-11, -1e-16])
    p.rezero(epsilon=1e-10)
    assert_not_equal(p.x, 0.0)
    assert_equal(p.y, 0.0)
    assert_equal(p.z, 0.0)


def test3():
    p = Point([1.0, 0.0, 0.0])
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([0, 1, 0])))
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([-1, 0, 0])))
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([0, -1, 0])))
    p.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(p, np.array([1, 0, 0])))
    p.rotate(np.pi/2, rot_axis=[0, 0, 1], rot_point=[1, 1, 0])
    assert_true(np.allclose(p, np.array([2, 1, 0])))


def test4():
    p0 = Point(np.zeros(3))
    p1 = Point(np.ones(3))
    p2 = Point(2 * np.ones(3))
    pm2 = Point(-2 * np.ones(3))
    p3 = Point([1, 2, 3])
    pm3 = Point([-1, -2, -3])
    plist = Points([p2, p1, p0, p3, pm2, pm3])
    plist = [p2, p1, p0, p3, pm2, pm3]
    plist.sort()
    assert_equal(plist.index(p0), 2)
    assert_true(p0 == plist[2])
    assert_equal(plist.index(p1), 3)
    assert_true(p1 == plist[3])
    assert_equal(plist.index(p2), 5)
    assert_true(p2 == plist[5])
    assert_equal(plist.index(pm2), 0)
    assert_true(pm2 == plist[0])
    assert_true(plist.index(p3), 4)
    assert_true(plist.index(pm3), 1)


def test5():
    p2 = Point(None, nd=2)
    p3 = Point(None, nd=3)
    assert_equal(p2, Point(np.zeros(2)))
    assert_equal(p3, Point(np.zeros(3)))


def test6():
    p2 = Point(np.zeros(3), nd=2)
    p3 = Point(np.zeros(2), nd=3)
    assert_equal(p2, Point(np.zeros(2)))
    assert_equal(p3, Point(np.zeros(3)))
    assert_true(len(p2) == 2)
    assert_true(len(p3) == 3)


def test7():
    p1 = Point([1, 0, 0])
    p2 = Point([0, 1, 0])
    p3 = Point([0, 0, 1])
    pts = Points([p1, p2, p3])
    assert_equal(pts[0], p1)
    assert_equal(pts[1], p2)
    assert_equal(pts[2], p3)
    [assert_equal(pts[i], pi) for i, pi in enumerate([p1, p2, p3])]
    [assert_is_instance(pt, Point) for pt in pts]

    assert_true(np.allclose(pts.x, [1, 0, 0]))
    assert_true(np.allclose(pts.y, [0, 1, 0]))
    assert_true(np.allclose(pts.z, [0, 0, 1]))

    pts[0] = np.ones(3)
    assert_true(np.allclose(pts[0], np.ones(3)))

    pts[1] = Point(np.ones(3))
    assert_true(np.allclose(pts[1], Point(np.ones(3))))

    assert_true(np.allclose(pts.x, [1, 1, 0]))
    assert_true(np.allclose(pts.y, [1, 1, 0]))
    assert_true(np.allclose(pts.z, [1, 1, 1]))

    pts.z = np.zeros(3)
    assert_true(np.allclose(pts.x, [1, 1, 0]))
    assert_true(np.allclose(pts.z, [0, 0, 0]))


def test8():
    p1 = Point([1, 0, 0])
    p2 = Point([0, 1, 0])
    p3 = Point([0, 0, 1])
    pts = Points([p1, p2, p3])

    p1.z = 1
    assert_equal(pts[0], p1)
    assert_true(pts[0] is p1)


def test9():
    v1 = Point(np.ones(3))
    v2 = Point(np.arange(3) + 1)
    v3 = Point(3 * np.ones(3))
    pts = Points([v1, v2, v3])
    assert_is_instance(pts, Points)
    [assert_is_instance(p2, Point) for p2 in pts]


def test10():
    p1 = Point([1, 0, 0])
    p2 = Point([0, 1, 0])
    p3 = Point([0, 0, 1])
    pts = Points([p1, p2, p3])
    assert_equal(pts[0], p1)
    assert_equal(pts[1], p2)
    assert_equal(len(pts[:3]), 3)
    [assert_equal(pts[i], pi) for i, pi in enumerate([p1, p2, p3])]


def test11():
    v0 = Point(np.arange(10))
    assert_equal(v0.argmax(), len(v0)-1)


def test12():
    pts = Points([[1, 1, 1], [1, 0, 1], [0, 1, 0]])
    pts += Point([1, 1, 1])
    assert_equal(pts[0], Point([2, 2, 2]))
    assert_equal(pts[1], Point([2, 1, 2]))
    assert_equal(pts[2], Point([1, 2, 1]))


def test13():
    pts1 = Points([[1, 1, 1], [1, 0, 1], [0, 1, 0]])
    pts2 = pts1 + Point([1, 1, 1])
    assert_equal(pts2[0], Point([2, 2, 2]))
    assert_equal(pts2[1], Point([2, 1, 2]))
    assert_equal(pts2[2], Point([1, 2, 1]))


def test14():
    pts1 = Points([[1, 1, 1], [1, 0, 1], [0, 1, 0]])
    pts2 = Point([1, 1, 1]) + pts1
    assert_equal(pts2[0], Point([2, 2, 2]))
    assert_equal(pts2[1], Point([2, 1, 2]))
    assert_equal(pts2[2], Point([1, 2, 1]))


def test15():
    pts1 = Points([[1, 1, 1], [1, 0, 1], [0, 1, 0]])
    pts2 = [1, 1, 1] + pts1
    assert_equal(pts2[0], Point([2, 2, 2]))
    assert_equal(pts2[1], Point([2, 1, 2]))
    assert_equal(pts2[2], Point([1, 2, 1]))


def test16():
    pts = Points([[1, 1, 1], [1, 0, 1], [0, 1, 0]])
    pts -= [1, 1, 1]
    assert_equal(pts[0], Point([0, 0, 0]))
    assert_equal(pts[1], Point([0, -1, 0]))
    assert_equal(pts[2], Point([-1, 0, -1]))


def test17():
    pts1 = Points([[1, 1, 1], [1, 0, 1], [0, 1, 0]])
    pts2 = Point([1, 1, 1]) - pts1
    assert_equal(pts2[0], Point([0, 0, 0]))
    assert_equal(pts2[1], Point([0, -1, 0]))
    assert_equal(pts2[2], Point([-1, 0, -1]))


def test18():
    pts1 = Points([[1, 1, 1], [1, 0, 1], [0, 1, 0]])
    pts2 = [1, 1, 1] - pts1
    assert_equal(pts2[0], Point([0, 0, 0]))
    assert_equal(pts2[1], Point([0, -1, 0]))
    assert_equal(pts2[2], Point([-1, 0, -1]))


def test19():
    p1 = Point([1, 0, 0])
    p2 = Point([0, 1, 0])
    p3 = Point([0, 0, 1])
    pts = Points([p1, p2, p3])

    pts += Point([1, 1, 1])
    [assert_equal(pts[i], pi) for i, pi in enumerate([p1, p2, p3])]


def test20():
    p1 = Point([1, 0, 0])
    p2 = Point([0, 1, 0])
    p3 = Point([0, 0, 1])
    pts = Points([p1, p2, p3])

    pts -= Point([1, 1, 1])
    [assert_equal(pts[i], pi) for i, pi in enumerate([p1, p2, p3])]


if __name__ == '__main__':
    nose.runmodule()
