#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.core.atoms import BasisAtom
from sknano.core.crystallography import Crystal2DLattice, Crystal3DLattice


def test1():
    from sknano.core.atoms import Atom
    xatom = BasisAtom()
    assert_is_instance(xatom, (Atom, BasisAtom))


def test2():
    lattice = Crystal2DLattice.square(a=1.0)
    elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
    for element in elements:
        xatom = BasisAtom(element, lattice=lattice)
        assert_equals(xatom.element, element)
        assert_equals(xatom.m, xatom.mass)

    xatom = BasisAtom()
    for c in ('x', 'y', 'z'):
        assert_equals(getattr(xatom, c), 0.0)


def test3():
    lattice = Crystal2DLattice.square(a=1.0)
    atom = BasisAtom('C', lattice=lattice)
    print(atom)


if __name__ == '__main__':
    nose.runmodule()
