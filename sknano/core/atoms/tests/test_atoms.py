#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import Atom, Atoms
from sknano.testing import generate_atoms
from six.moves import range


def test_instantiation():
    atoms = Atoms()
    assert_is_instance(atoms, Atoms)
    for Z in range(1, 101):
        atoms.append(Atom(Z))
    assert_equals(atoms.Natoms, 100)


def test_delete():
    atoms = generate_atoms(elements='periodic_table')
    print(atoms)


if __name__ == '__main__':
    nose.runmodule()
