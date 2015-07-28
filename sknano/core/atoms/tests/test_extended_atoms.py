#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.core.atoms import XAtom, XAtoms
# from sknano.testing import generate_atoms
# from sknano.core.geometric_regions import Ellipsoid


def test_list_methods():
    xatoms1 = XAtoms()
    for Z in range(100, 0, -1):
        xatoms1.append(XAtom(Z))
    xatoms1.sort(key=lambda a: a.Z)
    xatoms2 = XAtoms()
    for Z in range(1, 101):
        xatoms2.append(XAtom(Z))
    assert_equal(xatoms1, xatoms2)


def test_generator_atoms():
    from sknano.generators import SWNTGenerator
    n, m = 5, 5
    nz = 2
    swnt = SWNTGenerator(n=n, m=m, nz=nz)
    assert_equals(swnt.N, 2 * n)
    assert_equals(swnt.Natoms, 2 * swnt.N * swnt.nz)
    assert_equals(swnt.Natoms_per_unit_cell, 2 * swnt.N)
    assert_equals(swnt.nz, nz)
    assert_equals(swnt.Natoms_per_tube, swnt.Natoms_per_unit_cell * swnt.nz)
    atoms = swnt.atoms
    assert_equals(atoms.Natoms, swnt.Natoms_per_tube)


# def test_atom_selections():
#     atoms = \
#         generate_atoms(generator_class='SWNTGenerator', n=5, m=5, nz=2)
#     #assert_true(a200 is atoms[20])
#     #a200NN = atoms.select_within(Ellipsoid(center=a200.r, r=2.5))
#     #assert_equals(a200NN.Natoms, 4)
#     #atoms.select(


if __name__ == '__main__':
    nose.runmodule()
