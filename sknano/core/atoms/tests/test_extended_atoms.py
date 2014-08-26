#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import Atoms, XAtom, XAtoms


def test_atoms():
    xatoms = XAtoms()
    assert_is_instance(xatoms, (Atoms, XAtoms))
    for Z in range(1, 101):
        xatoms.append(XAtom(Z))

    assert_equals(xatoms.Natoms, 100)


if __name__ == '__main__':
    nose.runmodule()
