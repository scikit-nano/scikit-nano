#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import XAtom, XAtoms


def test_instantiation():
    from sknano.core.atoms import Atoms
    xatoms = XAtoms()
    assert_is_instance(xatoms, (Atoms, XAtoms))


def test_list_methods():
    xatoms = XAtoms()
    for Z in range(1, 101):
        xatoms.append(XAtom(Z))
    xatoms.sort()


def test_atom_tree():
    Catoms = XAtoms()
    for Z in range(1, 101):
        Catoms.append(XAtom(Z))
    assert_equals(Catoms.Natoms, 100)


def test_structure_analysis():
    from sknano.generators import SWNTGenerator
    n, m = 10, 10
    nz = 10
    swnt = SWNTGenerator(n=n, m=m, nz=nz)
    assert_equals(swnt.N, 2 * n)
    assert_equals(swnt.Natoms, 2 * swnt.N)
    assert_equals(swnt.nz, nz)
    assert_equals(swnt.Natoms_per_tube, swnt.Natoms * swnt.nz)
    atoms = swnt.atoms
    assert_equals(atoms.Natoms, swnt.Natoms_per_tube)


if __name__ == '__main__':
    nose.runmodule()
