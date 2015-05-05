#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import CrystalLattice, \
    CrystalStructure
from sknano.core.atoms import Atom, Atoms, XAtom, XAtoms
from sknano.core.math import Point


def test_crystal_lattice():
    lattice = CrystalLattice()
    assert_is_instance(lattice, CrystalLattice)


def test_crystal_structure():
    xtal = CrystalStructure(None)
    assert_is_instance(xtal, CrystalStructure)


def test3():
    xtal = CrystalLattice(a=4.0, b=8.0, c=2.0, alpha=90,
                          beta=90, gamma=120)
    print(xtal)


def test4():
    xtal = CrystalStructure(basis=XAtoms(atoms=[XAtom('Au'), XAtom('Au')]),
                            a=4.0, b=8.0, c=2.0, alpha=90,
                            beta=90, gamma=120)
    print(xtal)


def test5():
    xtal = CrystalStructure(basis=XAtom(element='Au'),
                            a=4.0, b=8.0, c=2.0, alpha=90,
                            beta=90, gamma=120)
    print(xtal)


def test6():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)
    p = [2.1, 0.9, 0.5]
    assert_true(np.allclose(latt.wrap_fractional_coordinate(p),
                Point((0.1, 0.9, 0.5))))


def test7():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)

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

    recip_latt = CrystalLattice(a1=latt.b1, a2=latt.b2, a3=latt.b3)

    assert_equal(latt.a1, recip_latt.b1)
    assert_equal(latt.a2, recip_latt.b2)
    assert_equal(latt.a3, recip_latt.b3)


if __name__ == '__main__':
    nose.runmodule()
