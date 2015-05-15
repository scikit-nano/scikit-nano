#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import range

import nose
from nose.tools import *
from sknano.core.atoms import BasisAtom, BasisAtoms
from sknano.testing import generate_atoms
#from sknano.utils.geometric_shapes import Ellipsoid


def test_instantiation():
    from sknano.core.atoms import Atoms
    basis_atoms = BasisAtoms()
    assert_is_instance(basis_atoms, (Atoms, BasisAtoms))
    swnt_atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=1)
    for atom in swnt_atoms:
        basis_atoms.append(atom)

    assert_equal(len(basis_atoms), len(swnt_atoms))

    atoms = BasisAtoms(atoms=basis_atoms)
    assert_equal(len(atoms), len(swnt_atoms))

    atoms = BasisAtoms(atoms=basis_atoms.data)
    assert_equal(len(atoms), len(swnt_atoms))

    atoms = BasisAtoms(atoms=atoms)
    assert_equal(len(atoms), len(swnt_atoms))


def test_list_methods():
    basis_atoms1 = BasisAtoms()
    for Z in range(100, 0, -1):
        basis_atoms1.append(BasisAtom(Z))
    basis_atoms1.sort(key=lambda a: a.Z)
    basis_atoms2 = BasisAtoms()
    for Z in range(1, 101):
        basis_atoms2.append(BasisAtom(Z))
    assert_equal(basis_atoms1, basis_atoms2)


def test_generator_atoms():
    from sknano.generators import SWNTGenerator
    n, m = 10, 10
    nz = 10
    swnt = SWNTGenerator(n=n, m=m, nz=nz)
    assert_equals(swnt.N, 2 * n)
    assert_equals(swnt.Natoms, 2 * swnt.N * swnt.nz)
    assert_equals(swnt.Natoms_per_unit_cell, 2 * swnt.N)
    assert_equals(swnt.nz, nz)
    assert_equals(swnt.Natoms_per_tube, swnt.Natoms_per_unit_cell * swnt.nz)
    atoms = swnt.atoms
    assert_equals(atoms.Natoms, swnt.Natoms_per_tube)

    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    assert_equals(atoms.Natoms, 400)

    a100 = atoms.get_atom(100)
    assert_true(a100 is atoms[99])
    assert_equals(atoms.index(a100), 99)
    a200 = atoms.get_atom(200)
    assert_true(a200 is atoms[199])
    assert_equals(atoms.index(a200), 199)
    a300 = atoms.get_atom(300)
    assert_true(a300 is atoms[299])
    assert_equals(atoms.index(a300), 299)


def test_atom_selections():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    #assert_true(a200 is atoms[20])
    #a200NN = atoms.select_within(Ellipsoid(center=a200.r, r=2.5))
    #assert_equals(a200NN.Natoms, 4)
    #atoms.select(


if __name__ == '__main__':
    nose.runmodule()
