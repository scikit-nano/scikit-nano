#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.core.atoms import Trajectory
#from sknano.testing import generate_atoms


def test_instantiation():
    traj = Trajectory()


if __name__ == '__main__':
    nose.runmodule()
