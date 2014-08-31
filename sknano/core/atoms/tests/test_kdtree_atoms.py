#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import KDTAtom, KDTAtoms
from sknano.testing import generate_atoms
#from sknano.utils.geometric_shapes import Ellipsoid


def test_instantiation():
    from sknano.core.atoms import Atoms, XAtoms
    atoms = KDTAtoms()
    assert_is_instance(atoms, (Atoms, XAtoms, KDTAtoms))
    swnt_atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=1)
    for atom in swnt_atoms:
        atoms.append(atom)

    assert_equal(len(atoms), len(swnt_atoms))

    atoms = KDTAtoms(atoms=atoms)
    assert_equal(len(atoms), len(swnt_atoms))

    atoms = KDTAtoms(atoms=atoms.data)
    assert_equal(len(atoms), len(swnt_atoms))

    atoms = KDTAtoms(atoms=atoms)
    assert_equal(len(atoms), len(swnt_atoms))


def test_list_methods():
    atoms1 = KDTAtoms()
    for Z in range(100, 0, -1):
        atoms1.append(KDTAtom(Z=Z))
    atoms1.sort(key=lambda a: a.Z)
    atoms2 = KDTAtoms()
    for Z in range(1, 101):
        atoms2.append(KDTAtom(Z=Z))
    assert_equal(atoms1, atoms2)


def test_atom_tree():
    Catoms = KDTAtoms()
    for Z in range(1, 101):
        Catoms.append(KDTAtom(Z))
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

    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    assert_equals(atoms.Natoms, 400)

    atoms = KDTAtoms(atoms=atoms)
    assert_equals(atoms.Natoms, 400)

    atoms.assign_unique_ids()
    atoms.kNN = 3
    atoms.NN_cutoff = 2.0
    NNatoms = atoms.nearest_neighbors
    assert_equals(len(NNatoms), atoms.Natoms)

    atomCN = atoms.coordination_numbers
    assert_equals(len(atomCN), atoms.Natoms)

    a100 = atoms.get_atom(atomID=100)
    assert_true(a100 is atoms[99])
    assert_equals(atoms.index(a100), 99)
    a200 = atoms.get_atom(atomID=200)
    assert_true(a200 is atoms[199])
    assert_equals(atoms.index(a200), 199)
    a300 = atoms.get_atom(atomID=300)
    assert_true(a300 is atoms[299])
    assert_equals(atoms.index(a300), 299)


def test_atom_selections():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    a200 = atoms.get_atom(atomID=200)
    assert_true(a200 is atoms[199])
    assert_true(a200 == atoms[199])
    #a200NN = atoms.select_within(Ellipsoid(center=a200.r, r=2.5))
    #assert_equals(a200NN.Natoms, 4)


if __name__ == '__main__':
    nose.runmodule()
