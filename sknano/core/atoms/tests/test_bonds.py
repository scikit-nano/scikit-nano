#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_true
from pkg_resources import resource_filename
import numpy as np
from sknano.generators import SWNTGenerator
from sknano.io import DATAReader
from sknano.core.math import Vector
from sknano.testing import AtomsTestFixture
from sknano.core.atoms import get_angle, Bonds
# from sknano.core import dedupe, flatten


def test1():
    infile = resource_filename('sknano', 'data/nanotubes/1005_5cells.data')
    data = DATAReader(infile)
    atoms = data.atoms
    atoms.assign_unique_ids()
    atoms.update_attrs()
    bonds = atoms.bonds
    assert_equal(bonds.Nbonds, atoms.coordination_numbers.sum())
    print(bonds.mean_length)
    print(bonds.atoms.Natoms)
    # print('bonds.mean_angle: {}'.format(bonds.mean_angle))


def test2():
    atoms = SWNTGenerator(n=10, m=10, nz=3).atoms
    atoms.assign_unique_ids()
    atoms.update_attrs()


def test3():
    atoms = SWNTGenerator(n=20, m=0, nz=2).atoms
    assert_equal(atoms.Natoms, 160)
    atoms.assign_unique_ids()
    atoms.update_attrs()
    print(np.degrees(atoms.bonds.mean_length))
    print(np.degrees(atoms.bonds.mean_angle))
    atom0 = atoms[0]
    atom0bonds = atom0.bonds
    print('atom0: {}'.format(atom0))
    print('atom0.r: {}'.format(atom0.r))
    for NN in atom0.NN:
        print('NN.r: {}'.format(NN.r))
    print(atom0bonds.atoms.Natoms)
    print(atom0bonds.atoms.CM)
    print('atoms.bonds.Nbonds: {}'.format(atoms.bonds.Nbonds))
    print('atoms.bonds.atoms.Natoms: {}'.format(atoms.bonds.atoms.Natoms))


def test4():
    atoms = SWNTGenerator(10, 5, nz=2).atoms
    atoms.assign_unique_ids()
    atoms.update_attrs()
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


class TestCase(AtomsTestFixture):
    def test1(self):
        atoms = self.atoms
        atoms.update_attrs()
        for n in range(1, 4):
            print('n={}'.format(n))
            bonds = atoms[11].bonds[:n]
            # print(bonds)
            angles = bonds.angles
            bond_angle_pairs = bonds.bond_angle_pairs
            if bonds is not None:
                print('Nbonds: {}'.format(bonds.Nbonds))
                assert_equal(bonds.Nbonds, n)
            if bond_angle_pairs is not None:
                print('angles: {}'.format(angles))
                print('bond_angle_pairs: {}'.format(bond_angle_pairs))
                for b1, b2 in bond_angle_pairs:
                    pair_bonds = Bonds(bonds=[b1, b2])
                    bond_atoms = pair_bonds.atoms
                    print(bond_atoms.ids)
                    print(bond_atoms)
                    assert_equal(bond_atoms.Natoms, 3)
                    # assert_equal(get_angle(bond_atoms), pair_bonds.angles)

            if n == 1:
                assert_true(angles is None)
                assert_true(bond_angle_pairs is None)
            elif n == 2:
                assert_true(isinstance(angles, float))
                assert_true(len(bond_angle_pairs) == 1)
            else:
                assert_true(isinstance(angles, np.ndarray))
                assert_true(len(bond_angle_pairs) == 3)

if __name__ == '__main__':
    nose.runmodule()
