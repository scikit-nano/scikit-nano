#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_is_instance
from sknano.core.molecules import Molecule, Molecules, \
    BasisMolecule, BasisMolecules, LatticeMolecule, LatticeMolecules, \
    PBCMolecule, PBCMolecules


def test_instantiation():

    for mi in (Molecule, BasisMolecule, LatticeMolecule, PBCMolecule):
        m = mi()
        print(m)
        assert_is_instance(m, Molecule)

    for mi in (Molecules, BasisMolecules, LatticeMolecules, PBCMolecules):
        m = mi()
        print(m)
        assert_is_instance(m, Molecules)


if __name__ == '__main__':
    nose.runmodule()
