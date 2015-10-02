#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# import numpy as np

import nose
from nose.tools import assert_equal, assert_true, assert_is_instance
from sknano.core.atoms import Atom, Atoms
from sknano.core.refdata import element_symbols
from sknano.testing import generate_atoms


def test1():
    atoms = Atoms(verbose=True)
    assert_is_instance(atoms, Atoms)
    for Z in range(1, 101):
        atoms.append(Atom(Z))
    assert_equal(atoms.Natoms, 100)
    assert_equal(Atoms(atoms=atoms).Natoms, atoms.Natoms)
    assert_equal(Atoms(atoms=atoms.data).Natoms, atoms.Natoms)


def test2():
    atoms = generate_atoms(elements='periodic_table')
    atoms.assign_unique_ids()

    atoms_slice = atoms[5:10]
    print(atoms_slice)


def test3():
    atoms = generate_atoms(elements='periodic_table')
    atoms.assign_unique_ids()
    print(atoms[:5])

    atoms_slice = atoms[5:10]

    atoms[:5] = atoms_slice
    assert_equal(atoms[:5], atoms[5:10])
    print(atoms[:5])

    atoms[:5] = ['C', 'H', 'N', 'Ar', 'He']
    print(atoms[:8])

    atoms[0] = 'Au'
    print(atoms[:2])


def test4():
    atoms = generate_atoms(elements='periodic_table')
    atoms.assign_unique_ids()
    a1 = atoms[:10]
    a2 = atoms[:5]
    assert_equal(a1 + a2, atoms[:10])

    a1 = atoms[:5]
    a2 = atoms[:10]
    assert_equal(a1 + a2, atoms[:10])

    assert_equal(atoms + atoms.__atom_class__('H', id=1), atoms)
    assert_equal(atoms.__atom_class__('H', id=1) + atoms, atoms)

    a1 = atoms[:25]
    a2 = atoms[25:]
    assert_equal((a1 + a2).elements.tolist(), element_symbols)
    assert_equal((a1 + a2), atoms)


def test5():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=5, m=5, nz=5,
                       verbose=True)
    atoms.assign_unique_ids()
    for i, atom in enumerate(atoms):
        assert_true(atoms[i] is atoms.get_atom(atom.id))
        assert_equal(i, atoms.index(atom))


def test6():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=3, nz=10)
    atoms.assign_unique_ids()
    atoms.center_centroid()
    atoms.kNN = 6
    atoms.NNrc = 9.0
    new_atoms = atoms.filtered((atoms.z >= -5) & (atoms.z <= 5))
    assert_equal(atoms.kNN, 6)
    assert_equal(atoms.NNrc, 9.0)
    assert_equal(atoms.kNN, new_atoms.kNN)
    assert_equal(atoms.NNrc, new_atoms.NNrc)


def test7():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=5, m=0, nz=5)
    atoms.update_attrs()
    filtered = atoms.filtered_ids(list(range(5, 26)))
    atoms.filter_ids(list(range(5, 26)))
    assert_equal(filtered, atoms)
    assert_true(atoms.Natoms, 20)


def test8():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=5, m=0, nz=5)
    atoms.update_attrs()
    filtered1 = atoms.filtered(atoms.coordination_numbers > 1)
    filtered2 = atoms.filtered("coordination_numbers > 1")
    atoms.filter(atoms.coordination_numbers > 1)
    assert_equal(filtered1, atoms)
    assert_equal(filtered1, filtered2)


def test_instantiation():
    a = Atom('C')
    assert_is_instance(a, Atom)


def test_attributes():
    atom = Atom('C')
    assert_equal(atom.element, 'C')
    assert_equal(atom.Z, 6)


def test_comparisons():
    B = Atom('B')
    C = Atom('C')
    N = Atom('N')
    assert_true(B < C)
    assert_true(N > C)
    assert_true(B < N)


def test_property_changes():
    CAtom = Atom('C')
    CAtom.Z = 54
    XeAtom = Atom('Xe')
    assert_equal(CAtom, XeAtom)


def test_str():
    atom = Atom('C')
    print(atom)

if __name__ == '__main__':
    nose.runmodule()
