#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import assert_equal, assert_equals, assert_false, \
    assert_true, assert_is_instance

import numpy as np

import sknano.core.atoms
from sknano.core.atoms import StructureAtom, StructureAtoms
from sknano.generators import SWNTGenerator
from sknano.io import DATAReader
from sknano.structures import compute_Natoms
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
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_attrs()
        print(atoms.ids)
        for atom in atoms:
            print('atom: {}, bond.lengths: {}'.format(
                  atom.id, atom.bonds.lengths))

    def test3(self):
        atoms = generate_atoms(elements='periodic_table')
        atoms.assign_unique_ids()
        atoms.kNN = 3
        atoms_cp = atoms.copy()
        assert_equal(atoms.kNN, atoms_cp.kNN)

    def test4(self):
        atoms = \
            generate_atoms(generator_class='SWNTGenerator', n=3, m=0, nz=5)
        assert_equal(compute_Natoms((3, 0), nz=5), atoms.Natoms)
        assert_equal(atoms.Natoms, atoms.ids[-1])

    def test5(self):
        atoms = self.atoms
        assert_true(np.allclose(atoms.coords, atoms.atom_tree.data))
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_attrs()
        assert_equals(len(atoms.nearest_neighbors), atoms.Natoms)
        assert_equals(len(atoms.coordination_numbers), atoms.Natoms)

        # atoms.kNN = 3
        # atoms.NNrc = 2.0
        # atoms.update_attrs()
        # print(atoms.ids)

    def test6(self):
        atoms = self.atoms
        atoms.update_attrs()
        assert_equal(atoms.filtered(atoms.coordination_numbers == 1).Natoms,
                     atoms.coordination_counts[1])
        assert_equal(atoms.filtered(atoms.coordination_numbers == 3).Natoms,
                     atoms.coordination_counts[3])

    def test7(self):
        atoms = self.atoms
        atoms.update_attrs()
        assert_true(np.allclose(atoms.coordination_numbers,
                                atoms.neighbor_counts(2.0)))

    def test8(self):
        atoms = self.atoms
        atoms.update_attrs()
        # print(atoms.bonds.lengths)
        # print(atoms.neighbor_distances)
        assert_true(np.allclose(atoms.bonds.lengths,
                                atoms.neighbor_distances))

    def test9(self):
        atoms = self.atoms
        atoms.update_attrs()
        assert_true(np.allclose(atoms.volume, atoms.bounds.volume))

    def test10(self):
        atom = StructureAtom(element='C')
        assert_equal(atom.CN, 0)
        atom.CN = 3
        assert_equal(atom.CN, 3)

    def test_list_methods(self):
        atoms1 = StructureAtoms()
        for Z in range(100, 0, -1):
            atoms1.append(StructureAtom(Z=Z))
        atoms1.sort(key=lambda a: a.Z)
        atoms2 = StructureAtoms()
        for Z in range(1, 101):
            atoms2.append(StructureAtom(Z=Z))
        assert_equal(atoms1, atoms2)

    def test_atom_bonds(self):
        atoms = \
            generate_atoms(generator_class='SWNTGenerator', n=10, m=5, nz=1)
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.compute_POAVs()
        bonds = atoms.bonds
        assert_equal(len(bonds), atoms.coordination_numbers.sum())
        assert_equal(bonds.Nbonds, atoms.coordination_numbers.sum())

        for i, atom in enumerate(atoms):
            if atom.bonds.Nbonds > 1:
                print('atom.bonds.angles:\n'
                      '{}'.format(np.degrees(atom.bonds.angles)))
                for j, bond in enumerate(atom.bonds):
                    assert_true(np.allclose(bond.vector,
                                            atom.bonds.vectors[j]))
                    assert_equal(bond.length, atom.bonds.lengths[j])

    def test_structure_data(self):
        fname = resource_filename('sknano', 'data/nanotubes/1005_5cells.data')
        swnt = SWNTGenerator(n=10, m=5, nz=5)
        swnt_atoms = swnt.atoms
        swnt_atoms.compute_POAVs()
        data = DATAReader(fname)
        atoms = data.atoms
        atoms.compute_POAVs()
        assert_equals(swnt_atoms.Natoms, atoms.Natoms)

    def test_POAVs(self):
        atoms = \
            generate_atoms(generator_class='SWNTGenerator', n=5, m=0, nz=5)
        atoms.compute_POAVs()
        atoms.filter(atoms.coordination_numbers == 3)
        atom = atoms[10]
        assert_equals(atom.bonds.Nbonds, 3)
        for POAV in ('POAV1', 'POAV2', 'POAVR'):
            atom_POAV = getattr(atom, POAV)
            assert_is_instance(atom_POAV, getattr(sknano.core.atoms, POAV))

            sigma_pi_angles = np.degrees(atom_POAV.sigma_pi_angles)
            assert_false(np.all(np.isclose(sigma_pi_angles, 3 * [np.nan],
                                           equal_nan=True)))
            pyramidalization_angles = \
                np.degrees(atom_POAV.pyramidalization_angles)
            assert_false(np.all(np.isclose(pyramidalization_angles,
                                           3 * [np.nan],
                                           equal_nan=True)))
            misalignment_angles = \
                np.degrees(atom_POAV.misalignment_angles)
            assert_false(np.all(np.isclose(misalignment_angles,
                                           3 * [np.nan],
                                           equal_nan=True)))
            print(getattr(atom, POAV))

    def test_POAV_angles(self):
        atoms = \
            generate_atoms(generator_class='SWNTGenerator', n=10, m=0, nz=2)
        # atoms.NNrc = 2.0
        atoms.assign_unique_ids()
        atoms.compute_POAVs()

        for i, atom in enumerate(atoms):
            print('atom{}: {}'.format(atom.id, atom))
            for POAV in ('POAV1', 'POAV2', 'POAVR'):
                if getattr(atom, POAV) is not None:
                    atom_POAV = getattr(atom, POAV)
                    sigma_pi_angles = np.degrees(atom_POAV.sigma_pi_angles)
                    assert_false(np.all(np.isclose(sigma_pi_angles,
                                                   3 * [np.nan],
                                                   equal_nan=True)))
                    print('atom{}.{}.sigma_pi_angles:\n{}'.format(
                        atom.id, POAV, sigma_pi_angles))

                    pyramidalization_angles = \
                        np.degrees(atom_POAV.pyramidalization_angles)
                    print('atom{}.{}.pyramidalization_angles:\n{}'.format(
                        atom.id, POAV, pyramidalization_angles))
                    assert_false(np.all(np.isclose(pyramidalization_angles,
                                                   3 * [np.nan],
                                                   equal_nan=True)))

                    misalignment_angles = \
                        np.degrees(atom_POAV.misalignment_angles)
                    print('atom{}.{}.misalignment_angles:\n{}\n'.format(
                        atom.id, POAV, misalignment_angles))
                    assert_false(np.all(np.isclose(misalignment_angles,
                                                   3 * [np.nan],
                                                   equal_nan=True)))


if __name__ == '__main__':
    nose.runmodule()
