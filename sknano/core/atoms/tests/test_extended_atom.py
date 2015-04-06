#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.core.atoms import XAtom


def test_instantiation():
    from sknano.core.atoms import Atom
    xatom = XAtom()
    assert_is_instance(xatom, (Atom, XAtom))


def test_attributes():
    elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
    for element in elements:
        xatom = XAtom(element=element)
        assert_equals(xatom.element, element)
        assert_equals(xatom.m, xatom.mass)

    xatom = XAtom()
    for c in ('x', 'y', 'z'):
        assert_equals(getattr(xatom, c), 0.0)


def test_str():
    atom = XAtom('C')
    print(atom)


if __name__ == '__main__':
    nose.runmodule()
