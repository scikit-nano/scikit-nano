#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
#from sknano.generators import SWNTGenerator
from sknano.testing import generate_atoms


def test1():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=3, m=3, nz=5)

    atoms.assign_unique_ids()
    atoms.kNN = 3
    atoms.NNrc = 2.0
    NNlist = atoms.nearest_neighbors
    CNlist = atoms.coordination_numbers
    for i, NN in enumerate(NNlist):
        assert_equals(NN.NNN, CNlist[i])


if __name__ == '__main__':
    nose.runmodule()
