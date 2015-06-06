#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import Crystal2DLattice, Reciprocal2DLattice
from sknano.core.atoms import Atom, Atoms, XAtom, XAtoms
from sknano.core.math import Point, transformation_matrix
from sknano.core.refdata import aCC, dVDW


def test1():
    latt = Crystal2DLattice(a=4.0, b=8.0, gamma=120)
    print(latt)
    assert_equal(str(latt),
                 "Crystal2DLattice(a=4.0, b=8.0, gamma=120.0)")


def test2():
    a = np.sqrt(3) * 1.42
    latt = Crystal2DLattice(a=a, b=a, gamma=120)
    hexlatt = Crystal2DLattice.hexagonal(a)
    assert_equal(latt, hexlatt)


def test3():
    a = np.sqrt(3) * 1.42
    latt = Crystal2DLattice(a=a, b=a, gamma=90)
    square = Crystal2DLattice.square(a)
    assert_equal(latt, square)


def test4():
    a = np.sqrt(3) * 1.42
    latt = Crystal2DLattice(a=a, b=a, gamma=60)
    a1 = latt.a1
    a2 = latt.a2

    rotated_a1 = a1.copy()
    rotated_a2 = a2.copy()
    xfrm = transformation_matrix(angle=-np.pi/6)
    rotated_a1.rotate(transform_matrix=xfrm)
    rotated_a2.rotate(transform_matrix=xfrm)

    latt.rotate(angle=-np.pi/6)

    assert_equal(latt.a1, rotated_a1)
    assert_equal(latt.a2, rotated_a2)
    assert_true(np.allclose(latt.orientation_matrix, xfrm))

    rotated_latt = Crystal2DLattice(a1=rotated_a1, a2=rotated_a2)
    assert_equal(rotated_a1, rotated_latt.a1)
    assert_equal(rotated_a2, rotated_latt.a2)
    assert_true(np.allclose(latt.orientation_matrix,
                            rotated_latt.orientation_matrix))


def test5():
    a = np.sqrt(3) * aCC
    latt = Crystal2DLattice(a=a, b=a, gamma=60)
    recip_latt = \
        Reciprocal2DLattice(a_star=latt.reciprocal_lattice.a_star,
                            b_star=latt.reciprocal_lattice.b_star,
                            gamma_star=latt.reciprocal_lattice.gamma_star)
    assert_equal(latt, recip_latt.reciprocal_lattice)
    assert_equal(latt.reciprocal_lattice, recip_latt)


def test6():
    a = np.sqrt(3) * aCC
    l1 = Crystal2DLattice.square(a)
    l2 = Crystal2DLattice.square(2*a)
    assert_true(l1 < l2)


if __name__ == '__main__':
    nose.runmodule()
