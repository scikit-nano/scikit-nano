#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
#from sknano.core.atoms import Atom
from sknano.core.molecules import Molecule


def test_instantiation():
    m = Molecule()
    assert_is_instance(m, Molecule)


if __name__ == '__main__':
    nose.runmodule()
