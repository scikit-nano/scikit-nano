#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import UnitCell, Crystal3DLattice, \
    Crystal3DStructure, AlphaQuartz, DiamondStructure, HexagonalStructure, \
    FCCStructure, Gold, Copper, CaesiumChlorideStructure, \
    RocksaltStructure, ZincblendeStructure, MoS2
from sknano.core.atoms import BasisAtoms


def test1():
    diamond = DiamondStructure()
    assert_true(isinstance(diamond.basis, BasisAtoms))
    assert_equal(diamond.basis.Natoms, 8)


def test2():
    gold = Gold()
    assert_true(isinstance(gold.basis, BasisAtoms))
    assert_equal(gold.basis.Natoms, 4)


def test3():
    copper = Copper()
    assert_true(isinstance(copper.basis, BasisAtoms))
    assert_equal(copper.basis.Natoms, 4)


def test4():
    quartz = AlphaQuartz()
    assert_true(isinstance(quartz.basis, BasisAtoms))
    assert_equal(quartz.basis.Natoms, 9)


def test5():
    caesium_chloride = CaesiumChlorideStructure()
    assert_true(isinstance(caesium_chloride.basis, BasisAtoms))
    assert_equal(caesium_chloride.basis.Natoms, 2)


def test6():
    rock_salt = RocksaltStructure()
    assert_true(isinstance(rock_salt.basis, BasisAtoms))
    assert_equal(rock_salt.basis.Natoms, 8)


def test7():
    zincblende = ZincblendeStructure()
    assert_true(isinstance(zincblende.basis, BasisAtoms))
    assert_equal(zincblende.basis.Natoms, 8)


def test8():
    molydisulphide = MoS2()
    print(molydisulphide)
    assert_true(isinstance(molydisulphide.basis, BasisAtoms))
    assert_equal(molydisulphide.basis.Natoms, 6)

if __name__ == '__main__':
    nose.runmodule()
