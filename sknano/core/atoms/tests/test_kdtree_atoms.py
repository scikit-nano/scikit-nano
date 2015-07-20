#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
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


def test10():
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


if __name__ == '__main__':
    nose.runmodule()
