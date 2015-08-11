#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equals
from sknano.testing import AtomsTestFixture, generate_atoms


class TestCase(AtomsTestFixture):

    def test1(self):
        atoms = self.atoms

        atoms.kNN = 6
        atoms.NNrc = 9.0
        for atom in atoms:
            assert_equals(atom.CN, 0)
            atom.CN = 3
            assert_equals(atom.CN, 3)

        atoms.update_attrs()
        atoms = atoms.filtered((atoms.z >= -5) & (atoms.z <= 5))
        print('Natoms: {}'.format(atoms.Natoms))
        for atom in atoms:
            assert_equals(atom.CN, atoms.kNN)

    def test2(self):
        atoms = self.atoms
        atoms.update_attrs()
        print(atoms.neighbors.tolist())
        print(atoms.neighbor_distances)

    def test3(self):
        atoms = \
            generate_atoms(generator_class='SWNTGenerator', n=10, m=0, nz=5)
        atoms.kNN = 30
        atoms.NNrc = 10
        atoms.update_neighbors(cutoffs=[1.5, 2.5, 2.9, 3.8, 4.3])
        # print(atoms.neighbors.tolist())
        print(atoms.neighbor_distances.tolist())
        # print('first_neighbors:\n{}'.format(atoms.first_neighbors))
        # print('second_neighbors:\n{}'.format(atoms.second_neighbors))
        print('len(first_neighbors): {}'.format(
              [len(atom.first_neighbors) for atom in atoms]))
        print('len(second_neighbors): {}'.format(
              [len(atom.second_neighbors) for atom in atoms]))
        print('len(third_neighbors): {}'.format(
              [len(atom.third_neighbors) for atom in atoms]))
        print('len(4th_neighbors): {}'.format(
              [len(atom.get_neighbors(4)) for atom in atoms]))
        print('len(5th_neighbors): {}'.format(
              [len(atom.get_neighbors(5)) for atom in atoms]))


if __name__ == '__main__':
    nose.runmodule()
