#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import KDTAtom


def test_instantiation():
    from sknano.core.atoms import Atom, XAtom
    kdtatom = KDTAtom()
    assert_is_instance(kdtatom, (Atom, XAtom, KDTAtom))


def test_attributes():
    elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
    for element in elements:
        kdtatom = KDTAtom(element=element)
        assert_equals(kdtatom.element, element)


if __name__ == '__main__':
    nose.runmodule()
