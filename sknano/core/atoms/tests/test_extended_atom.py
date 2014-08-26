#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import Atom, XAtom


def test_instantiation():
    xatom = XAtom()
    assert_is_instance(xatom, (Atom, XAtom))


if __name__ == '__main__':
    nose.runmodule()
