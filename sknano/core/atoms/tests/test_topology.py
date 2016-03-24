#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_true, assert_raises
import numpy as np
# from sknano.generators import SWNTGenerator
# from sknano.io import DATAReader
from sknano.core.math import Vector
from sknano.core.atoms import compute_angle, compute_dihedral, compute_improper
from sknano.testing import AtomsTestFixture


class Tests(AtomsTestFixture):

    def test1(self):
        atoms = self.atoms
        print('atoms.Natoms: {}'.format(atoms.Natoms))
        # atoms.set_pbc('z')
        # atoms.update_attrs()

        all_bonds = atoms.all_bonds
        assert_equal(all_bonds.Nbonds, atoms.coordination_numbers.sum())

        bonds = atoms.bonds
        assert_equal(all_bonds.Nbonds, 2 * bonds.Nbonds)

        neighbor_ids = atoms.get_bonds(neighbors_only=True, ids_only=True)
        print('len(neighbor_ids): {}'.format(len(neighbor_ids)))
        print('len(all_bonds.atom_ids): {}'.format(len(all_bonds.atom_ids)))
        print('len(bonds.atom_ids): {}'.format(len(bonds.atom_ids)))

    def test2(self):
        atoms = self.atoms
        all_angles = atoms.all_angles
        angles = atoms.angles
        assert_equal(all_angles.Nangles, len(set(all_angles.atom_ids)))
        print('\nall_angles.Nangles: {}'.format(all_angles.Nangles))
        print('angles.Nangles: {}'.format(angles.Nangles))

    def test3(self):
        atoms = self.atoms
        all_dihedrals = atoms.all_dihedrals
        dihedrals = atoms.dihedrals
        assert_equal(all_dihedrals.Ndihedrals,
                     len(set(all_dihedrals.atom_ids)))
        print('\nall_dihedrals.Ndihedrals: {}'.format(all_dihedrals.Ndihedrals))
        print('dihedrals.Ndihedrals: {}'.format(dihedrals.Ndihedrals))

        # for angle in dihedrals.angles:
        #     print(np.degrees(angle))

    def test4(self):
        atoms = self.atoms
        all_impropers = atoms.all_impropers
        impropers = atoms.impropers
        assert_equal(all_impropers.Nimpropers,
                     len(set(all_impropers.atom_ids)))
        print('\aall_impropers.Nimpropers: {}'.format(all_impropers.Nimpropers))
        print('impropers.Nimpropers: {}'.format(impropers.Nimpropers))

    def test5(self):
        atoms = self.atoms
        atoms.set_pbc('xyz')
        assert_equal(atoms.Natoms, 100)

        # print(np.degrees(atoms.bonds.mean_length))
        # print(np.degrees(atoms.bonds.mean_angle))
        atom0 = atoms[0]
        # atom0bonds = atom0.bonds
        # print('atom0: {}'.format(atom0))
        print('atom0.r: {}'.format(atom0.r))
        for NN in atom0.NN:
            print('NN.r: {}'.format(NN.r))
        assert_equal(atom0.bonds.Nbonds, 1)
        # assert_equal(atom0.bonds.atoms.Natoms, 2)
        # print(atom0bonds.atoms.CM)
        # print('atoms.bonds.Nbonds: {}'.format(atoms.bonds.Nbonds))

    def test6(self):
        atoms = self.atoms
        atoms.update_neighbors()
        atoms.center_centroid()
        print(atoms.centroid)
        bond = atoms.get_atom(atoms.Natoms // 2).bonds[0]
        print(bond)
        print('bond.atoms.coords:\n{}'.format(bond.atoms.coords))
        bond_centroid = bond.centroid
        print('bond.centroid: {}'.format(bond.centroid))
        rot_axis = Vector(p0=[0, 0, bond_centroid.z], p=bond_centroid.p)
        bond.rotate(angle=np.pi/2, axis=rot_axis)
        print(bond)
        print('bond.atoms.coords:\n{}'.format(bond.atoms.coords))
        assert_true(np.allclose(bond_centroid, bond.centroid))

    def test7(self):
        atoms = self.atoms
        atoms.update_neighbors()
        index = 11
        for n in range(1, 4):
            print('n={}'.format(n))
            # print(atoms[index].bonds)
            bonds = atoms[index].bonds[:n]
            # print(bonds)
            # angles = atoms.angles
            # bond_angle_pairs = bonds.bond_angle_pairs
            if bonds is not None:
                print('Nbonds: {}'.format(bonds.Nbonds))
                assert_equal(bonds.Nbonds, n)
            # if bond_angle_pairs is not None:
            #     print('angles: {} degrees'.format(np.degrees(angles)))
            #     # print('bond_angle_pairs: {}'.format(bond_angle_pairs))
            #     for b1, b2 in bond_angle_pairs:
            #         assert_equal(b1.origin, b2.origin)
            #         # print(b1.atom_ids)
            #         # print(b2.atom_ids)
            #         pair_bonds = Bonds(bonds=[b1, b2])

            # if n == 1:
            #     assert_true(angles is None)
            #     assert_true(bond_angle_pairs is None)
            # elif n == 2:
            #     assert_true(isinstance(angles, float))
            #     assert_true(len(bond_angle_pairs) == 1)
            # else:
            #     assert_true(isinstance(angles, np.ndarray))
            #     assert_true(len(bond_angle_pairs) == 3)

    def test8(self):
        atoms = self.atoms
        print(atoms.topology_stats)
        atoms.angles_in_degrees = True
        print(atoms.topology_stats)

if __name__ == '__main__':
    nose.runmodule()
