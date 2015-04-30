#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import CrystalLattice, \
    CrystalStructure
from sknano.core.atoms import Atom, Atoms, XAtom, XAtoms


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



if __name__ == '__main__':
    nose.runmodule()
