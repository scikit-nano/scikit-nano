#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename

import numpy as np

from sknano.core.atoms import Trajectory
from sknano.io import DUMPReader
#from sknano.testing import generate_atoms


def test1():
    traj = Trajectory()
    assert_equal(traj.Nsnaps, 0)


def test2():
    dump = \
        DUMPReader(resource_filename('sknano',
                                     'data/lammpstrj/0500_29cells.dump'))
    dump.remap_atomattr_names({'c_peratom_pe': 'pe', 'c_peratom_ke': 'ke'})

    traj = dump.trajectory
    traj.reference_snapshot = traj[0]
    prev_ss_atom = None
    for ss in traj:
        atoms = ss.atoms
        atoms = atoms.filter((atoms.z <= 20) & (atoms.z >= -20))
        atoms.update_attrs()
        atoms = atoms.filter((atoms.z <= 15) & (atoms.z >= -15))
        atom = atoms.get_atom(259)
        print('atom.NN:\n{}'.format(atom.NN))
        print('atom.bonds:\n{}\n'.format(atom.bonds))

        if prev_ss_atom is not None:
            assert_true(np.all(atom.NN.ids == prev_ss_atom.NN.ids))

        prev_ss_atom = atom


if __name__ == '__main__':
    nose.runmodule()
