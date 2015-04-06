#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

#import numpy as np

from sknano.core.refdata import element_symbols
from sknano.testing import generate_atoms


def test_generate_atoms():
    atoms = generate_atoms(elements=element_symbols)
    assert_equals(atoms.Natoms, len(element_symbols))
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    assert_equals(atoms.Natoms, 400)

if __name__ == '__main__':
    nose.runmodule()
