#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import numpy as np

from sknano.core.crystallography import Crystal2DLattice, Crystal2DStructure
from sknano.core.atoms import BasisAtoms


def test1():
    lattice = Crystal2DLattice(a=3, b=3, gamma=60)
    lattice.rotate(angle=-np.pi/6)
    structure = Crystal2DStructure(lattice=lattice, basis=['C', 'C'],
                                   coords=[[0, 0, 0], [1.5, 0, 0]],
                                   cartesian=True)
    assert_true(isinstance(structure.basis, BasisAtoms))
    assert_equal(structure.basis.Natoms, 2)


if __name__ == '__main__':
    nose.runmodule()
