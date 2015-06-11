#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import UnitCell, Crystal2DLattice, \
    Crystal2DStructure, Crystal3DLattice, Crystal3DStructure, \
    DiamondStructure, HexagonalStructure, FCCStructure, \
    Gold, Copper
from sknano.core.atoms import BasisAtom, BasisAtoms


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
