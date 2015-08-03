#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import CrystalLattice, ReciprocalLattice
# from sknano.core.atoms import Atom, Atoms, XAtom, XAtoms
from sknano.core.math import Point, transformation_matrix
from sknano.core.refdata import aCC, element_data


def test1():
    latt = CrystalLattice(a=4.0, b=8.0, c=2.0, alpha=90,
                          beta=90, gamma=120)
    assert_equal(latt.a, 4.0)
    assert_equal(latt.b, 8.0)
    assert_equal(latt.c, 2.0)
    assert_equal(latt.alpha, 90)
    assert_equal(latt.beta, 90)
    assert_equal(latt.gamma, 120)


def test2():
    a = np.sqrt(3) * aCC
    latt = CrystalLattice(a=a, b=a, c=2*element_data['C']['VanDerWaalsRadius'],
                          alpha=90, beta=90, gamma=120)
    print(latt)
    a1 = latt.a1
    a2 = latt.a2
    a3 = latt.a3

    xfrm = transformation_matrix(angle=-np.pi/6)

    rotangle = -np.pi / 6
    for v in (a1, a2, a3):
        v.rotate(angle=rotangle)

    latt.rotate(angle=rotangle, axis='z')
    print(latt)

    assert_equal(latt.a1, a1)
    assert_equal(latt.a2, a2)
    assert_equal(latt.a3, a3)

    assert_true(np.allclose(latt.orientation_matrix, xfrm))


def test3():
    a = np.sqrt(3) * aCC
    latt = CrystalLattice(a=a, b=a, c=2*element_data['C']['VanDerWaalsRadius'],
                          alpha=90, beta=90, gamma=120)
    print(latt)
    hex_latt = \
        CrystalLattice.hexagonal(a, 2 * element_data['C']['VanDerWaalsRadius'])
    print(hex_latt)
    assert_equal(latt, hex_latt)


def test4():
    a = np.sqrt(3) * aCC
    latt = CrystalLattice(a=a, b=a, c=2*element_data['C']['VanDerWaalsRadius'],
                          alpha=90, beta=90, gamma=120)
    print(latt)
    recip_latt = \
        ReciprocalLattice(a_star=latt.reciprocal_lattice.a_star,
                          b_star=latt.reciprocal_lattice.b_star,
                          c_star=latt.reciprocal_lattice.c_star,
                          alpha_star=latt.reciprocal_lattice.alpha_star,
                          beta_star=latt.reciprocal_lattice.beta_star,
                          gamma_star=latt.reciprocal_lattice.gamma_star)
    assert_equal(latt, recip_latt.reciprocal_lattice)
    assert_equal(latt.reciprocal_lattice, recip_latt)


def test5():
    a = np.sqrt(3) * aCC
    cubic_latt = CrystalLattice(a=a, b=a, c=a, alpha=90, beta=90, gamma=90)
    assert_equal(cubic_latt, CrystalLattice.cubic(a))


def test6():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)
    print(latt)
    p = [2.1, 0.9, 0.5]
    assert_true(np.allclose(latt.wrap_fractional_coordinate(p),
                Point((0.1, 0.9, 0.5))))


def test7():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)
    print(latt)

    a = latt.a1
    b = latt.a2
    c = latt.a3
    G = np.matrix([[a.dot(a), a.dot(b), a.dot(c)],
                   [b.dot(a), b.dot(b), b.dot(c)],
                   [c.dot(a), c.dot(b), c.dot(c)]])

    assert_true(np.allclose(latt.metric_tensor, G))


def test8():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)
    print(latt)

    recip_latt = CrystalLattice(a1=latt.b1, a2=latt.b2, a3=latt.b3)
    print(recip_latt)

    assert_equal(latt.a1, recip_latt.b1)
    assert_equal(latt.a2, recip_latt.b2)
    assert_equal(latt.a3, recip_latt.b3)


def test9():
    a = np.sqrt(3) * aCC
    assert_true(CrystalLattice.cubic(a) < CrystalLattice.cubic(2*a))


if __name__ == '__main__':
    nose.runmodule()
