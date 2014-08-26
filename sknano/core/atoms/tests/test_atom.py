#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import Atom


def test_atom_attributes():
    a = Atom('C')
    assert_is_instance(a, Atom)
    assert_equals(a.element, 'C')

    for c in ('x', 'y', 'z'):
        assert_equals(getattr(a, c), 0.0)


if __name__ == '__main__':
    nose.runmodule()
