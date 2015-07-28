#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

import numpy as np

from sknano.structures import compute_Natoms
from sknano.testing import generate_atoms


def test1():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=5, m=0, nz=5)
    atoms.assign_unique_ids()
    atoms.update_attrs()
    assert_true(np.allclose(atoms.volume, atoms.bounds.volume))


if __name__ == '__main__':
    nose.runmodule()
