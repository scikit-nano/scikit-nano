#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

import numpy as np

from sknano.structures import compute_Natoms
from sknano.testing import generate_atoms


def test_NN_parameters():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=3, nz=10)
    atoms.assign_unique_ids()
    atoms.center_CM()
    atoms.kNN = 6
    atoms.NNrc = 9.0
    new_atoms = atoms.filtered((atoms.z >= -5) & (atoms.z <= 5))
    assert_equal(atoms.kNN, 6)
    assert_equal(atoms.NNrc, 9.0)
    assert_equal(atoms.kNN, new_atoms.kNN)
    assert_equal(atoms.NNrc, new_atoms.NNrc)


def test_filtered_ids():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=5, m=0, nz=5)
    atoms.update_attrs()
    filtered = atoms.filtered_ids(list(range(5, 26)))
    atoms.filter_ids(list(range(5, 26)))
    assert_equal(filtered, atoms)
    assert_true(atoms.Natoms, 20)


def test_filtered():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=5, m=0, nz=5)
    atoms.update_attrs()
    filtered1 = atoms.filtered(atoms.coordination_numbers > 1)
    filtered2 = atoms.filtered("coordination_numbers > 1")
    atoms.filter(atoms.coordination_numbers > 1)
    assert_equal(filtered1, atoms)
    assert_equal(filtered1, filtered2)


def test1():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=3, nz=10)
    atoms.assign_unique_ids()
    atoms.kNN = 6
    atoms.NNrc = 9.0
    for atom in atoms:
        assert_equals(atom.CN, 0)
        atom.CN = 3
        assert_equals(atom.CN, 3)

    atoms.update_attrs()
    atoms = atoms.filtered((atoms.z >= -5) & (atoms.z <= 5))
    print('Natoms: {}'.format(atoms.Natoms))
    for atom in atoms:
        assert_equals(atom.CN, atoms.kNN)


def test2():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=3, nz=10)
    atoms.assign_unique_ids()
    atoms.kNN = 3
    atoms.NNrc = 2.0
    atoms.update_attrs()
    print(atoms.ids)
    for atom in atoms:
        print('atom: {}, bond.lengths: {}'.format(atom.id, atom.bonds.lengths))


def test3():
    atoms = generate_atoms(elements='periodic_table')
    atoms.kNN = 3
    atoms_cp = atoms.copy()
    assert_equal(atoms.kNN, atoms_cp.kNN)


def test4():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=0, nz=5)
    assert_equal(compute_Natoms((3, 0), nz=5), atoms.Natoms)
    assert_equal(atoms.Natoms, atoms.ids[-1])


def test5():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=0, nz=5)
    assert_true(np.allclose(atoms.coords, atoms.atom_tree.data))

    # atoms.kNN = 3
    # atoms.NNrc = 2.0
    # atoms.update_attrs()
    # print(atoms.ids)


def test6():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=0, nz=5)
    atoms.update_attrs()
    assert_equal(atoms.filtered(atoms.coordination_numbers == 1).Natoms,
                 atoms.coordination_counts[1])
    assert_equal(atoms.filtered(atoms.coordination_numbers == 3).Natoms,
                 atoms.coordination_counts[3])


def test7():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=0, nz=5)
    atoms.update_attrs()
    assert_true(np.allclose(atoms.coordination_numbers,
                            atoms.neighbor_counts(2.0)))


if __name__ == '__main__':
    nose.runmodule()
