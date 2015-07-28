#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
import sknano.core.atoms
from sknano.core.atoms import POAVAtom
from sknano.testing import generate_atoms


def test_attributes():
    elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
    for element in elements:
        atom = POAVAtom(element=element)
        assert_equals(atom.element, element)


def test_POAVs():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=5, m=0, nz=5)
    atoms.compute_POAVs()
    atoms.filter(atoms.coordination_numbers == 3)
    atom = atoms[10]
    assert_equals(atom.bonds.Nbonds, 3)
    for POAV in ('POAV1', 'POAV2', 'POAVR'):
        assert_is_instance(getattr(atom, POAV),
                           getattr(sknano.core.atoms, POAV))
        print(getattr(atom, POAV))


if __name__ == '__main__':
    nose.runmodule()
