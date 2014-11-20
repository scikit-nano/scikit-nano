#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import POAVAtom, POAV1, POAV2, POAVR
from sknano.testing import generate_atoms


def test_instantiation():
    from sknano.core.atoms import Atom, XAtom, KDTAtom
    atom = POAVAtom()
    assert_is_instance(atom, (Atom, XAtom, KDTAtom, POAVAtom))


def test_attributes():
    elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
    for element in elements:
        atom = POAVAtom(element=element)
        assert_equals(atom.element, element)


def test_POAV1():
    #atom = POAVAtom(element='C')
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    atoms.update_nearest_neighbors()
    atoms.update_coordination_numbers()
    atoms.update_bonds()
    atom100 = atoms.get_atom(100)
    assert_equals(atom100.atomID, 100)
    assert_true(atom100.POAV1 is None)
    assert_equals(atom100.bonds.Nbonds, 3)
    setattr(atom100, 'POAV1', POAV1(atom100.bonds))


def test_POAV2():
    #atom = POAVAtom(element='C')
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    atoms.update_nearest_neighbors()
    atoms.update_coordination_numbers()
    atoms.update_bonds()
    atom100 = atoms.get_atom(100)
    assert_equals(atom100.atomID, 100)
    assert_true(atom100.POAV2 is None)
    assert_equals(atom100.bonds.Nbonds, 3)
    setattr(atom100, 'POAV2', POAV2(atom100.bonds))


def test_POAVR():
    #atom = POAVAtom(element='C')
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    atoms.update_nearest_neighbors()
    atoms.update_coordination_numbers()
    atoms.update_bonds()
    atom100 = atoms.get_atom(100)
    assert_equals(atom100.atomID, 100)
    assert_true(atom100.POAVR is None)
    assert_equals(atom100.bonds.Nbonds, 3)
    setattr(atom100, 'POAVR', POAVR(atom100.bonds))


if __name__ == '__main__':
    nose.runmodule()
