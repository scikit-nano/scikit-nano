#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import warnings
warnings.simplefilter('always')

import nose
from nose.tools import assert_equal, assert_true, assert_raises, assert_false
import numpy as np

from sknano.core.crystallography import Crystal2DLattice
from sknano.core.structures import Crystal2DStructure, \
    AlphaQuartz, BetaQuartz, DiamondStructure, HexagonalStructure, \
    BCCStructure, FCCStructure, Gold, Copper, CaesiumChlorideStructure, \
    RocksaltStructure, ZincblendeStructure, MoS2, Fullerene
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
    assert_equal(structure.Natoms, 4)
    structure2 = FCCStructure('Au', scaling_matrix=2)
    assert_equal(2 ** 3 * structure.Natoms, structure2.Natoms)


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


def test15():
    buckyball = Fullerene(60)
    fcc_buckyball = FCCStructure(structure=buckyball,
                                 coords=[[0, 0, 0]])
    assert_equal(buckyball.Natoms * 4, fcc_buckyball.Natoms)
    # print(fcc_buckyball)
    # print(fcc_buckyball.Natoms)
    # buckyball.center_centroid()
    # print(buckyball.centroid)
    # basis = buckyball.basis
    # coords = basis.coords
    # print(coords)
    # basis.center_centroid()
    # print(basis.centroid)

    # fcc_structure = FCCStructure(basis='C', lattice=buckyball.lattice)
    # print(fcc_structure.basis)
    # print(fcc_structure.basis.coords)
    # new_basis = BasisAtoms()

    # for i, tvec in enumerate(fcc_structure.basis[:].coords):
    #     # print(tvec)
    #     basis = buckyball.basis.copy()
    #     # print(basis)
    #     basis.translate(tvec)
    #     # print(basis.centroid)
    #     # fcc_structure.basis[i] = basis
    #     # print(fcc_structure.basis[i].Natoms)
    #     new_basis.extend(basis)
    # # print(fcc_structure.basis.Natoms)
    # # print(fcc_structure.basis[0])
    # # print(fcc_structure.basis)
    # # assert_false(fcc_structure.basis[0] == fcc_structure.basis[1])
    # print(new_basis.Natoms)


# def test16():
#     structure = BetaQuartz()
#     print(structure)
#     assert_true(isinstance(structure.basis, BasisAtoms))
#     assert_equal(structure.basis.Natoms, 9)


if __name__ == '__main__':
    nose.runmodule()
