#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import warnings

import nose
from nose.tools import *

import numpy as np

from sknano.core.math import e1, e2, e3, xhat, yhat, zhat, \
    Point, Vector, vector as vec, Rx, Ry, Rz, NullVector


def test1():
    v = Vector()
    assert_true(np.allclose(v, np.zeros(3)))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.zeros(3)))

    v = Vector([1, 1, 1])
    assert_true(np.allclose(v, np.ones(3)))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.ones(3)))

    v = Vector(p=[1, 1, 1])
    assert_true(np.allclose(v, np.ones(3)))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.ones(3)))

    v = Vector(p0=[1, 1, 1])
    assert_true(np.allclose(v, -np.ones(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, np.zeros(3)))

    v = Vector(p=[1, 1, 1], p0=[1, 1, 1])
    assert_true(np.allclose(v, np.zeros(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, np.ones(3)))


def test2():
    v = Vector()
    v.p0 = np.ones(3)
    assert_true(np.allclose(v, -np.ones(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, np.zeros(3)))

    v = Vector(p=[1., 1., 1.], p0=[1., 1., 1.])
    assert_true(np.allclose(v, np.zeros(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, np.ones(3)))
    v.p0 = np.zeros(3)
    assert_true(np.allclose(v, np.ones(3)))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.ones(3)))
    v.x = 5.0
    assert_equal(v.x, 5.0)
    assert_equal(v.x, v.p.x)
    v.y = -5.0
    assert_equal(v.y, -5.0)
    assert_equal(v.y, v.p.y)
    v.z = 0.0
    assert_equal(v.z, 0.0)
    assert_equal(v.z, v.p.z)

    v.p0 = np.array([0.5, -10.0, 2.5])
    assert_equal(v.p0.x, 0.5)
    assert_equal(v.p0.y, -10.0)
    assert_equal(v.p0.z, 2.5)


def test3():
    v = Vector(np.zeros(3))
    assert_true(np.allclose(v.p0, np.zeros(3)))

    v = Vector(np.ones(3))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.ones(3)))

    v = Vector(np.ones(3), p0=[1, 1, 1])
    assert_true(np.allclose(v, np.ones(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, 2*np.ones(3)))

    v[:] += np.ones(3)
    assert_true(np.allclose(v, 2 * np.ones(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, 3 * np.ones(3)))

    v = Vector()
    v[:] += np.ones(3)
    assert_true(np.allclose(v, np.ones(3)))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.ones(3)))

    v1 = Vector([1.0, 0.0], p0=[5., 5.])
    v2 = Vector([1.0, 1.0], p0=[5., 5.])
    v3 = v1 + v2
    #print('v3: {}'.format(v3))
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.array([2.0, 1.0])))
    assert_true(np.allclose(v3.p0, np.array([5.0, 5.0])))
    assert_true(np.allclose(v3.p, np.array([7.0, 6.0])))
    v3 = v2 + v1
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.array([2.0, 1.0])))
    assert_true(np.allclose(v3.p0, np.array([5.0, 5.0])))
    assert_true(np.allclose(v3.p, np.array([7.0, 6.0])))


def test4():
    v1 = Vector([1.0, 0.0])
    v2 = Vector([1.0, 1.0])

    c = np.dot(np.asarray(v1), np.asarray(v2))
    assert_equal(c, 1.0)

    c = vec.dot(v1, v2)
    assert_equal(c, 1.0)


def test5():
    v1 = Vector([1.0, 0.0])
    v2 = Vector([1.0, 1.0])

    v3 = vec.cross(v1, v2)
    assert_equal(v3, np.cross(np.asarray(v1), np.asarray(v2)))

    v1 = Vector([1, 0], p0=[1, 1])
    v2 = Vector([0, 1], p0=[1, 1])

    v3 = vec.cross(v1, v2)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))

    v1 = Vector([1, 0], p0=[1, 0])
    v2 = Vector([0, 1], p0=[0, 1])

    v3 = vec.cross(v1, v2)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))


def test6():
    v1 = Vector([1.0, 0.0, 0.0])
    v2 = Vector([1.0, 1.0, 0.0])

    c = np.dot(np.asarray(v1), np.asarray(v2))
    assert_equal(c, 1.0)

    c = vec.dot(v1, v2)
    assert_equal(c, 1.0)


def test7():
    v1 = Vector([1, 0, 0])
    v2 = Vector([0, 1, 0])
    v3 = vec.cross(v1, v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_true(np.allclose(v3.p0, v2.p0))

    v1 = Vector([1, 0, 0], p0=[1, 1, 1])
    v2 = Vector([0, 1, 0], p0=[1, 1, 1])

    v3 = vec.cross(v1, v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_true(np.allclose(v3.p0, v2.p0))

    v1 = Vector([1, 0, 0], p0=[1, 0, 0])
    v2 = Vector([0, 1, 0], p0=[0, 1, 0])

    v3 = vec.cross(v1, v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_false(np.allclose(v3.p0, v2.p0))

    v1 = Vector([1, 2, 3], p0=[1, 0, 0])
    v2 = Vector([4, 5, 6], p0=[0, 1, 0])
    v3 = vec.cross(v1, v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_false(np.allclose(v3.p0, v2.p0))


def test8():
    assert_true(np.allclose(vec.cross(e1, e2), e3))
    assert_equal(vec.cross(e1, e2), e3)

    assert_true(np.allclose(vec.cross(e1, e3), -e2))
    assert_equal(vec.cross(e1, e3), -e2)

    assert_true(np.allclose(vec.cross(e2, e3), e1))
    assert_equal(vec.cross(e2, e3), e1)


def test9():
    v = Vector([1.0, 0.0, 0.0])
    v.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(v, np.array([0, 1, 0])))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.array([0, 1, 0])))
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)
    v.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(v, np.array([-1, 0, 0])))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.array([-1, 0, 0])))
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)
    v.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(v, np.array([0, -1, 0])))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.array([0, -1, 0])))
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)
    v.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(v, np.array([1, 0, 0])))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.array([1, 0, 0])))
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)

    v = Vector(p0=[1, 1, 1])
    assert_true(np.allclose(v, -np.ones(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, np.array([0, 0, 0])))
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)
    v.rotate(np.pi/2, rot_axis='z')
    assert_true(np.allclose(v, np.array([1, -1, -1])))
    assert_true(np.allclose(v.p0, np.array([-1, 1, 1])))
    assert_true(np.allclose(v.p, np.zeros(3)))

    v = Vector(p0=[1, 1, 1])
    assert_true(np.allclose(v, -np.ones(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, np.array([0, 0, 0])))
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)
    v.rotate(np.pi/2, rot_axis='z', rot_point=[1, 1, 1])
    assert_true(np.allclose(v, np.array([1, -1, -1])))
    assert_true(np.allclose(v.p0, np.array([1, 1, 1])))
    assert_true(np.allclose(v.p, np.array([2, 0, 0])))
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)


def test10():
    v = Vector(p0=[1, 1, 1])
    assert_true(np.allclose(v, -np.ones(3)))
    assert_true(np.allclose(v.p0, np.ones(3)))
    assert_true(np.allclose(v.p, np.array([0, 0, 0])))
    v.rotate(np.pi/2, rot_axis='z', rot_point=[1, 1, 1])
    assert_is_instance(v.p0, Point)
    assert_is_instance(v.p, Point)
    assert_true(np.allclose(v, np.array([1, -1, -1])))
    assert_true(np.allclose(v.p0, np.array([1, 1, 1])))
    assert_true(np.allclose(v.p, np.array([2, 0, 0])))


def test11():
    v = Vector()
    v += np.ones(3)
    assert_true(np.allclose(v, np.ones(3)))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, np.ones(3)))

    v.p += np.ones(3)
    assert_true(np.allclose(v, 2 * np.ones(3)))
    assert_true(np.allclose(v.p0, np.zeros(3)))
    assert_true(np.allclose(v.p, 2 * np.ones(3)))

    v1 = Vector()
    dr = Vector(np.ones(3))
    v2 = v1 + dr
    #print('v2: {}'.format(v2))
    #print('v2.p: {}'.format(v2.p))
    assert_true(np.allclose(v2, np.ones(3)))
    assert_true(np.allclose(v2.p0, np.zeros(3)))
    assert_true(np.allclose(v2.p, np.ones(3)))


def test12():

    v1 = Vector()
    dr = Vector(np.ones(3))
    v3 = Vector()
    v3 = Vector(v1 + dr)
    #v3[:] = v1 + dr
    #print('v3: {}'.format(v3))
    #print('v3.p: {}'.format(v3.p))
    assert_true(np.allclose(v3, np.ones(3)))
    assert_true(np.allclose(v3.p0, np.zeros(3)))
    assert_true(np.allclose(v3.p, np.ones(3)))


def test13():

    u = Vector([5, 6, 7])
    assert_true(np.allclose(
        u.projection(Vector([1, 0, 0])), Vector([5, 0, 0])))
    assert_true(np.allclose(
        u.projection(Vector([1, 1, 1])), Vector([6, 6, 6])))


def test14():
    v = Vector([5, 5, 0])
    assert_almost_equal(v.unit_vector.length, 1.0)


def test15():
    u = Vector([1, 0])
    v = Vector([1, 1])
    assert_almost_equal(vec.angle(u, v), np.pi/4)


def test16():
    u = Vector([1, 0, 0])
    v = Vector([1, 1, 0])
    w = Vector([0, 1, 1])
    assert_equals(vec.scalar_triple_product(u, v, w), 1)


def test17():
    u = Vector([1, 0, 0])
    v = Vector([1, 1, 0])
    w = Vector([0, 1, 1])
    assert_is_instance(vec.vector_triple_product(u, v, w), Vector)


def test18():
    u = Vector([1, 2, 3])
    v = Vector([1, 2, 3])
    w = Vector([1, 1, 1])
    assert_equals(u, v)
    assert_true(u == v)
    assert_not_equals(u, w)
    assert_false(u == w)


def test19():
    u = Vector([10, 10, 17])
    v = Vector([10, 10, 17])
    assert_true(u == v)
    v.normalize()
    assert_true(np.allclose(u.unit_vector, v))


def test20():
    a = Vector([1, 1, 1])
    b = Vector([1, 0, 0])
    assert_equals(vec.scalar_projection(a, b), 1.0)

    a = Vector([1, 0, 0])
    b = Vector([1, 1, 1])
    assert_equals(vec.scalar_projection(a, b), 1.0 / np.sqrt(3))


def test21():
    a = Vector([1, 1, 1])
    b = Vector([1, 0, 0])
    assert_equals(vec.vector_projection(a, b), Vector([1, 0, 0]))

    a = Vector([1, 0, 0])
    b = Vector([1, 1, 1])
    assert_equals(vec.vector_projection(a, b), Vector(1/3 * np.ones(3)))

    a = Vector([5, 6, 7])
    b = Vector([1, 1, 1])
    assert_equals(vec.vector_projection(a, b), Vector(3 * [6]))


def test22():
    a = Vector([1, 1, 1])
    b = Vector([1, 0, 0])
    assert_equals(vec.vector_rejection(a, b), Vector([0, 1, 1]))

    a = Vector([1, 0, 0])
    b = Vector([1, 1, 1])
    assert_equals(vec.vector_rejection(a, b), Vector([2/3, -1/3, -1/3]))

    a = Vector([5, 6, 7])
    b = Vector([1, 1, 1])
    assert_equals(vec.vector_rejection(a, b), a - Vector(3 * [6]))


def test23():
    for v1, v2 in zip([e1, e2, e3], [xhat, yhat, zhat]):
        assert_true(np.allclose(v1, v2))
        assert_equals(v1, v2)


def test24():
    assert_true(np.allclose(Vector([5,0,0]), 5 * e1))


def test25():
    assert_equal(e1 * e1, e1.dot(e1))
    assert_equal(e1 * e1, 1.0)

    with assert_raises(ValueError):
        e1.dot(Vector([1.0, 0.0]))

    with assert_raises(ValueError):
        Vector([1.0, 0.0]).dot(e1)


def test26():
    v = Vector([5, 0, 0])
    assert_equal(5 * e1, v)
    assert_equal(e1 * 5, v)

    e1copy = e1.copy()
    assert_true(e1 == e1copy)
    assert_true(not e1 is e1copy)
    assert_true(not (e1 is e1copy))
    assert_true(e1 is not e1copy)
    e1copy *= 5
    assert_equal(5 * e1, e1copy)
    assert_equal(e1copy, v)


def test27():
    v = Vector([5, 0, 0])
    v += 5
    assert_equal(v, Vector([10, 5, 5]))

    v += v
    assert_equal(v, Vector([20, 10, 10]))
    #with assert_raises(Val


def test28():
    v = Vector([5, 0, 0])
    v -= 5
    assert_equal(v, Vector([0, -5, -5]))


def test29():
    u = Vector([5, 0, 0])
    v = u.copy()

    u *= v

    assert_equal(u, 25)


def test30():
    u = Vector([5, 0, 0])
    v = u.copy()

    with assert_raises(TypeError):
        v /= v

    assert_equal(u, v)

    #v /= np.array([2.0])
    #assert_equal(u / 2.0, v)


def test31():
    u = Vector([5, 0, 0])
    v = u.copy()

    with assert_raises(TypeError):
        v //= v

    assert_equal(u, v)

    #v //= np.array([2.0])
    #assert_equal(u // 2.0, v)


def test32():
    u = Vector([5, 0, 0])
    v = u.copy()

    with assert_raises(TypeError):
        v **= v

    assert_equal(u, v)

    #v **= np.array([2.0])
    #assert_equal(u**2, v)


def test33():
    v1 = Vector([1, 0, 0])
    v2 = Vector([0, 1, 0])
    v3 = v1.cross(v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_true(np.allclose(v3.p0, v2.p0))

    v1 = Vector([1, 0, 0], p0=[1, 1, 1])
    v2 = Vector([0, 1, 0], p0=[1, 1, 1])

    v3 = v1.cross(v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_true(np.allclose(v3.p0, v2.p0))

    v1 = Vector([1, 0, 0], p0=[1, 0, 0])
    v2 = Vector([0, 1, 0], p0=[0, 1, 0])

    v3 = v1.cross(v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_false(np.allclose(v3.p0, v2.p0))

    v1 = Vector([1, 2, 3], p0=[1, 0, 0])
    v2 = Vector([4, 5, 6], p0=[0, 1, 0])
    v3 = v1.cross(v2)
    assert_is_instance(v3, Vector)
    assert_true(np.allclose(v3, np.cross(np.asarray(v1), np.asarray(v2))))
    assert_true(np.allclose(v3.p0, v1.p0))
    assert_false(np.allclose(v3.p0, v2.p0))


def test34():
    assert_equal(e1.angle(e1), 0.0)
    assert_equal(e1.angle(e2), np.pi/2)
    assert_equal(e1.angle(-e1), np.pi)


def test35():
    v0 = Vector(np.arange(3) + 2)
    I = np.identity(3)
    v1 = Vector((v0.row_matrix * I).A.flatten(), p0=v0.p0)
    v2 = Vector((I * v0.column_matrix).A.flatten(), p0=v0.p0)

    assert_true(np.allclose(v0, v1))
    assert_true(np.allclose(v0, v2))


def test36():
    v0 = Vector(np.arange(10))
    assert_equal(v0.argmax(), len(v0)-1)


def test37():
    assert_equal(Vector([0, 0, 0]), NullVector())
    a = NullVector()
    b = NullVector()
    assert_equal(a, b)
    assert_true(a is b)


if __name__ == '__main__':
    nose.runmodule()
