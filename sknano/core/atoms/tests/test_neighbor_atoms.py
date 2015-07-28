#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.testing import generate_atoms


def test1():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=3, m=0, nz=5)
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
