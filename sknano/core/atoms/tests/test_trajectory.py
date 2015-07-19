#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from pkg_resources import resource_filename

import numpy as np

from sknano.core.atoms import Trajectory
from sknano.io import DUMPReader
# from sknano.testing import generate_atoms


def test1():
    traj = Trajectory()
    assert_equal(traj.Nsnaps, 0)


def test2():
    dump = \
        DUMPReader(resource_filename('sknano',
                                     'data/lammpstrj/0500_29cells.dump'),
                   attrmap={'c_peratom_pe': 'pe', 'c_peratom_ke': 'ke'})

    traj = dump.trajectory
    traj.reference_snapshot = traj[0]
    prev_ss_atom = None
    for i, ss in enumerate(traj, start=1):
        atoms = ss.atoms
        atoms = atoms.filtered((atoms.z <= 20) & (atoms.z >= -20))
        atoms.update_attrs()
        atoms = atoms.filtered((atoms.z <= 15) & (atoms.z >= -15))
        atom = atoms.get_atom(259)
        # print('ss={}, atom.NN:\n{}'.format(i, atom.NN))
        # print('ss={}, atom.bonds:\n{}\n'.format(i, atom.bonds))

        if prev_ss_atom is not None:
            print('\ntimestep: {}'.format(i))
            # print('reference_atom')
            # print(prev_ss_atom.reference_atom)
            # print(atom.reference_atom)
            # print(prev_ss_atom.NN)
            # print(atom.NN)
            # print(prev_ss_atom.bonds)
            # print(atom.bonds)
            print(prev_ss_atom.bonds.lengths)
            print(atom.bonds.lengths)
            print(prev_ss_atom.NN.ids)
            print(atom.NN.ids)
            assert_true(np.all(atom.NN.ids == prev_ss_atom.NN.ids))

        prev_ss_atom = atom


if __name__ == '__main__':
    nose.runmodule()
