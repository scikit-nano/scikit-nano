#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename
from sknano.core.atoms import Trajectory
from sknano.io import DUMPData
#from sknano.testing import generate_atoms


def test1():
    traj = Trajectory()
    assert_equal(traj.Nsnaps, 0)


def test2():
    dumpdata = \
        DUMPData(resource_filename('sknano',
                                   'data/lammpstrj/0500_29cells.dump'))
    trajectory = dumpdata.trajectory
    trajectory.reference_snapshot = trajectory[0]
    ss2975_atoms = trajectory[0].atoms
    assert_equal(ss2975_atoms.Natoms, trajectory.reference_atoms.Natoms)


if __name__ == '__main__':
    nose.runmodule()
