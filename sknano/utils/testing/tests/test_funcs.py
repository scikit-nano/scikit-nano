#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *

import numpy as np

from sknano.utils.testing import generate_atoms


def test_generate_atoms():
    atoms = generate_atoms(from='periodic_table', filter=None)
    


if __name__ == '__main__':
    nose.runmodule()
