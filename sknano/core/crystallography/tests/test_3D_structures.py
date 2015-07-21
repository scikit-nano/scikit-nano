#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import UnitCell, Crystal3DLattice, \
    Crystal3DStructure, AlphaQuartz, DiamondStructure, HexagonalStructure, \
    BCCStructure, FCCStructure, Gold, Copper, CaesiumChlorideStructure, \
    RocksaltStructure, ZincblendeStructure, MoS2
from sknano.core.atoms import BasisAtoms


def test1():
    structure = DiamondStructure()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 8)


def test2():
    structure = Gold()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 4)


def test3():
    structure = Copper()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 4)


def test4():
    structure = AlphaQuartz()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 9)


def test5():
    structure = CaesiumChlorideStructure()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 2)


def test6():
    structure = RocksaltStructure()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 8)


def test7():
    structure = ZincblendeStructure()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 8)


def test8():
    structure = MoS2()
    print(structure)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 6)


def test9():
    with assert_raises(ValueError):
        BCCStructure()


def test10():
    s1 = FCCStructure('Au')
    print(s1)
    s2 = Gold()
    print(s2)
    assert_equal(s1.scaling_matrix, s2.scaling_matrix)
    assert_equal(s1.lattice, s2.lattice)
    assert_equal(s1.unit_cell, s2.unit_cell)
    assert_equal(s2, s1)


def test11():
    assert_equal(FCCStructure('Au'), FCCStructure(basis='Au'))


if __name__ == '__main__':
    nose.runmodule()
