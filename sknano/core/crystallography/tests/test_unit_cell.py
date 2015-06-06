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


if __name__ == '__main__':
    nose.runmodule()
