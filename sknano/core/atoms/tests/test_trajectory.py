#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename
import unittest

import nose
from nose.tools import assert_equal, assert_true

import numpy as np

from sknano.core.atoms import Trajectory
from sknano.io import DUMPReader
# from sknano.testing import generate_atoms


def test1():
    traj = Trajectory()
    assert_equal(traj.Nsnaps, 0)


class TrajectoryTestFixture(unittest.TestCase):

    def setUp(self):
        self.dump = \
            DUMPReader(resource_filename('sknano',
                                         'data/lammpstrj/0500_29cells.dump'),
                       attrmap={'c_peratom_pe': 'pe', 'c_peratom_ke': 'ke'})


class Tests(TrajectoryTestFixture):

    def test1(self):
        traj = self.dump.trajectory
        assert_true('{!s}'.format(traj).startswith('Trajectory(snapshots='))
        assert_true('{!s}'.format(traj[0]).startswith('Snapshot(trajectory='))

    def test2(self):
        traj = self.dump.trajectory
        assert_equal(traj.Nsnaps, self.dump.Nsnaps)

    def test3(self):
        traj = self.dump.trajectory
        assert_equal(traj.nselect, traj.nselected)
        assert_true(traj.aselect is traj.atom_selection)
        assert_true(traj.tselect is traj.time_selection)
        assert_true(traj.Nsnaps == len(traj.timesteps))
        traj.time_selection.one(2975)
        assert_true(len(traj.timesteps) == 1)

    def test4(self):
        traj = self.dump.trajectory
        traj.time_selection.one(2975)
        assert_true(len(traj.timesteps) == 1)
        self.dump.time_selection.all()
        assert_equal(traj.nselected, self.dump.nselected)
        assert_equal(self.dump.nselected, self.dump.Nsnaps)
        assert_equal(self.dump.Nsnaps, len(self.dump.timesteps))

    def test5(self):
        traj = self.dump.trajectory
        prev_ss_atom = None
        for ss in traj[10:-1:5]:
            atoms = ss.atoms
            # atoms = atoms.filtered((atoms.z <= 20) & (atoms.z >= -20))
            atoms = atoms.select("(z <= 20) and (z >= -20)")
            atoms.update_attrs()
            # atoms = atoms.filtered((atoms.z <= 15) & (atoms.z >= -15))
            atoms = atoms.select("(z <= 15) and (z >= -15)")
            # assert_equal(atoms, selection)
            atom = atoms.get_atom(259)
            # assert_equal(atom, selection.get_atom(259))
            # print('ss={}, atom.NN:\n{}'.format(i, atom.NN))
            # print('ss={}, atom.bonds:\n{}\n'.format(i, atom.bonds))

            if prev_ss_atom is not None:
                print('\ntimestep: {}'.format(ss.timestep))
                # print('reference_atom')
                # print(prev_ss_atom.reference_atom)
                # print(atom.reference_atom)
                # print(prev_ss_atom.NN)
                # print(atom.NN)
                # print(prev_ss_atom.bonds)
                # print(atom.bonds)
                print(prev_ss_atom.reference_atom_bonds.lengths)
                print(atom.reference_atom_bonds.lengths)
                print(prev_ss_atom.reference_atom_neighbors.ids)
                print(atom.reference_atom_neighbors.ids)
                assert_true(np.all(atom.reference_atom_neighbors.ids ==
                                   prev_ss_atom.reference_atom_neighbors.ids))

            prev_ss_atom = atom


if __name__ == '__main__':
    nose.runmodule()
