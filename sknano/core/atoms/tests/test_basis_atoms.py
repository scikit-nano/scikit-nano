#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import range

import nose
from nose.tools import *
from sknano.core.atoms import BasisAtom, BasisAtoms
from sknano.core.crystallography import Crystal2DLattice, Crystal3DLattice
from sknano.testing import generate_atoms
#from sknano.utils.geometric_shapes import Ellipsoid


def test1():
    lattice = Crystal3DLattice.cubic(a=5.0)
    basis = BasisAtoms(atoms=['C', 'C'], lattice=lattice)
    print(basis)


if __name__ == '__main__':
    nose.runmodule()
