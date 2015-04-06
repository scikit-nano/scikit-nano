#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
#from sknano.core.atoms import Atom, Atoms
from sknano.core.molecules import Molecules
#from sknano.testing import generate_atoms


def test_instantiation():
    molecules = Molecules()
    assert_is_instance(molecules, Molecules)


if __name__ == '__main__':
    nose.runmodule()
