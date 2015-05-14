#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import warnings

import nose
from nose.tools import *

import numpy as np

from sknano.core.math import Vector, Quaternion


def test1():
    with assert_raises(TypeError):
        H = Quaternion()

    with assert_raises(ValueError):
        H = Quaternion(np.ones(3))

    H = Quaternion(np.zeros(4))
    print(H)


def test2():
    H = Quaternion(np.arange(4))
    assert_equal(H.w, H.real)
    assert_equal(H.w, 0.)
    assert_is_instance(H.imag, list)
    assert_is_instance(H.v, Vector)
    assert_true(np.allclose(H.v, Vector(np.arange(3) + 1)))


def test3():
    q0 = Quaternion(np.arange(4) + 1)
    q1 = Quaternion(np.arange(4) + 2)
    q2 = Quaternion(np.arange(4) + 3)
    s = 7.0

    assert_true(np.abs(q0.norm - 5.477225575) < 1e-8)
    assert_true(np.abs(q1.norm - 7.348469228) < 1e-8)
    assert_true(np.abs(q2.norm - 9.273618495) < 1e-8)

    assert_equal(-q0, Quaternion(-(np.arange(4) + 1)))
    assert_equal(q0.conjugate,
                 Quaternion.from_components(w=q0.w, x=-q0.x, y=-q0.y, z=-q0.z))


if __name__ == '__main__':
    nose.runmodule()
