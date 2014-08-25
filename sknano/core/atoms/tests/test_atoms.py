#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
import unittest

from sknano.core.atoms import Atom, Atoms


class TestAtoms(unittest.TestCase):

    def test_atom(self):
        C_atom = Atom('C')
        print(C_atom)
        C_atom.x = 1.5
        print(C_atom)
        C_atom.rezero_coords()

    def test_atoms(self):
        atoms = Atoms()
        atom1 = Atom('C')
        atoms.append(atom1)
        print(atoms)
        atom2 = Atom('B')
        atoms.append(atom2)
        print(atoms)
        atom3 = Atom('N')
        atoms.append(atom3)
        print(atoms)
        atom3.x = 1.5
        print(atoms)


if __name__ == '__main__':
    unittest.main()
