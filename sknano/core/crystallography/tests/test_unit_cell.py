#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import UnitCell, Crystal2DLattice, \
    Crystal2DStructure, Crystal3DLattice, Crystal3DStructure
from sknano.core.atoms import BasisAtom, BasisAtoms
from sknano.core.math import Point
from sknano.core.refdata import aCC


def test1():
    a = np.sqrt(3) * aCC
    lattice = Crystal2DLattice(a=a, b=a, gamma=60)
    lattice.rotate(angle=-np.pi/6)
    basis = ['C', 'C']
    coords = [[0, 0, 0], [aCC, 0, 0]]
    structure = Crystal2DStructure(lattice, basis, coords, cartesian=True)
    print(structure.unit_cell)
    cell = UnitCell(lattice, basis, coords, cartesian=True)
    print(cell)


def test2():
    lattice = Crystal3DLattice.hexagonal(a=5.0, c=10.0)
    cell = UnitCell(lattice=lattice, basis=['C', 'C'],
                    coords=[[0, 0, 0], [1/2, 1/2, 1/2]])
    print(cell)
    print(cell.basis)


def test3():
    a = 5.0
    angle = 90
    basis = ['C', 'C']
    coords = [[0, 0, 0], [1/2, 1/2, 1/2]]
    latt1 = Crystal3DLattice.cubic(a)
    latt2 = \
        Crystal3DLattice(a=a, b=a, c=a, alpha=angle, beta=angle, gamma=angle)
    assert_equal(latt1, latt2)
    cell1 = UnitCell(lattice=latt1, basis=basis[:], coords=coords[:])
    cell2 = UnitCell(lattice=latt2, basis=basis[:], coords=coords[:])
    assert_equal(cell1, cell2)
    latt3 = Crystal3DLattice.cubic(10.0)
    assert_not_equal(latt1, latt3)


if __name__ == '__main__':
    nose.runmodule()
