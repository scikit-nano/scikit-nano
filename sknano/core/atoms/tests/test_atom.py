#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import Atom


def test_instantiation():
    a = Atom('C')
    assert_is_instance(a, Atom)


def test_attributes():
    atom = Atom('C')
    assert_equals(atom.element, 'C')
    for c in ('x', 'y', 'z'):
        assert_equals(getattr(atom, c), 0.0)
    assert_equals(atom.Z, 6)


def test_comparisons():
    B = Atom('B')
    C = Atom('C')
    N = Atom('N')
    assert_true(B < C)
    assert_true(N > C)
    assert_true(B < N)


if __name__ == '__main__':
    nose.runmodule()
