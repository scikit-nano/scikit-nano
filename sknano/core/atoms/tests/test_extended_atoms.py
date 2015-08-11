#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal
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
    assert_equal(swnt.N, 2 * n)
    assert_equal(swnt.Natoms, 2 * swnt.N * swnt.nz)
    assert_equal(swnt.Natoms_per_unit_cell, 2 * swnt.N)
    assert_equal(swnt.nz, nz)
    assert_equal(swnt.Natoms_per_tube, swnt.Natoms_per_unit_cell * swnt.nz)
    atoms = swnt.atoms
    assert_equal(atoms.Natoms, swnt.Natoms_per_tube)


# def test_atom_selections():
#     atoms = \
#         generate_atoms(generator_class='SWNTGenerator', n=5, m=5, nz=2)
#     #assert_true(a200 is atoms[20])
#     #a200NN = atoms.select_within(Ellipsoid(center=a200.r, r=2.5))
#     #assert_equal(a200NN.Natoms, 4)
#     #atoms.select(

def test_attributes():
    elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
    for element in elements:
        xatom = XAtom(element=element)
        assert_equal(xatom.element, element)
        assert_equal(xatom.m, xatom.mass)

    xatom = XAtom()
    for c in ('x', 'y', 'z'):
        assert_equal(getattr(xatom, c), 0.0)


if __name__ == '__main__':
    nose.runmodule()
