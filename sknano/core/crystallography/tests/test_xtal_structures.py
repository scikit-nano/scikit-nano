#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_true, assert_raises
import numpy as np

from sknano.core.crystallography import Crystal2DLattice, Crystal2DStructure, \
    AlphaQuartz, DiamondStructure, HexagonalStructure, \
    BCCStructure, FCCStructure, Gold, Copper, CaesiumChlorideStructure, \
    RocksaltStructure, ZincblendeStructure, MoS2
from sknano.core.atoms import BasisAtoms
from sknano.testing import generate_atoms, generate_structure


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
    assert_true(np.allclose(s1.scaling_matrix, s2.scaling_matrix))
    assert_equal(s1.lattice, s2.lattice)
    assert_equal(s1.unit_cell, s2.unit_cell)
    assert_equal(s2, s1)


def test11():
    assert_equal(FCCStructure('Au'), FCCStructure(basis='Au'))


def test12():
    structure = FCCStructure('Au')
    print(structure)
    structure2 = FCCStructure('Au', scaling_matrix=2)
    print(structure2)
    assert_equal(4 ** 3 * structure.Natoms, structure2.Natoms)


def test13():
    lattice = Crystal2DLattice(a=3, b=3, gamma=60)
    lattice.rotate(angle=-np.pi/6)
    structure = Crystal2DStructure(lattice=lattice, basis=['C', 'C'],
                                   coords=[[0, 0, 0], [1.5, 0, 0]],
                                   cartesian=True)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 2)


def test14():
    structure = \
        generate_structure(generator_class='SWNTGenerator', n=5, m=0, nz=2)
    print(structure)
    print(structure.crystal_cell)


if __name__ == '__main__':
    nose.runmodule()
