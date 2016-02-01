#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

# from collections import Counter
from operator import attrgetter

# import networkx as nx
import numpy as np

import nose
from nose.tools import assert_equal, assert_true, assert_false, \
    assert_is_instance

from sknano.core import flatten
from sknano.core.atoms import Atom, Atoms, BasisAtom, BasisAtoms, \
    StructureAtom, StructureAtoms
from sknano.core.crystallography import Crystal2DLattice, Crystal3DLattice
from sknano.core.math import Vectors
from sknano.core.refdata import element_symbols
from sknano.io import DATAReader
from sknano.testing import AtomsTestFixture


class Tests(AtomsTestFixture):

    def test1(self):
        a = Atom('C')
        assert_is_instance(a, Atom)

    def test2(self):
        atom = Atom('C')
        assert_equal(atom.element, 'C')
        assert_equal(atom.Z, 6)

        for element in element_symbols:
            atom = Atom(element=element)
            assert_equal(atom.element, element)

    def test3(self):
        B = Atom('B')
        C = Atom('C')
        N = Atom('N')
        assert_true(B < C)
        assert_true(C < N)
        assert_true(B < N)

    def test4(self):
        CAtom = Atom('C')
        CAtom.Z = 54
        XeAtom = Atom('Xe')
        assert_equal(CAtom, XeAtom)

    def test5(self):
        atoms1 = Atoms()
        for Z in range(100, 0, -1):
            atoms1.append(Atom(Z))
        atoms1.sort(key=lambda a: a.Z)

        atoms2 = Atoms()
        for Z in range(1, 101):
            atoms2.append(Atom(Z))
        assert_equal(atoms1, atoms2)

    def test6(self):
        atoms = self.periodic_table
        atoms_slice = atoms[5:10]
        assert_equal(len(atoms_slice), 5)
        # print(atoms_slice)

    def test7(self):
        atoms = self.periodic_table
        atoms_slice = atoms[5:10]
        atoms[:5] = atoms_slice
        assert_equal(atoms[:5], atoms[5:10])
        new_elements = ['C', 'H', 'N', 'Ar', 'He']
        atoms[:5] = new_elements
        [assert_equal(atom.element, element) for atom, element
         in zip(atoms[:5], new_elements)]
        atoms[0] = 'Au'
        assert_equal(atoms[0].element, 'Au')

    def test8(self):
        atoms = self.periodic_table
        a1 = atoms[:10]
        a2 = atoms[:5]
        assert_equal(a1 + a2, atoms[:10])

        a1 = atoms[:5]
        a2 = atoms[:10]
        assert_equal(a1 + a2, atoms[:10])

        assert_equal(atoms + atoms.__atom_class__('H', id=1), atoms)
        assert_equal(atoms.__atom_class__('H', id=1) + atoms, atoms)

        a1 = atoms[:25]
        a2 = atoms[25:]
        assert_equal((a1 + a2).elements.tolist(), element_symbols)
        assert_equal((a1 + a2), atoms)

    def test9(self):
        atoms = self.atoms
        for i, atom in enumerate(atoms):
            assert_true(atoms[i] is atoms.get_atom(atom.id))
            assert_equal(i, atoms.index(atom))

    def test10(self):
        atoms = Atoms(verbose=True)
        assert_is_instance(atoms, Atoms)
        for Z in range(1, 101):
            atoms.append(Atom(Z))
        assert_equal(atoms.Natoms, 100)
        assert_equal(Atoms(atoms=atoms).Natoms, atoms.Natoms)
        assert_equal(Atoms(atoms=atoms.data).Natoms, atoms.Natoms)

    def test11(self):
        atom = BasisAtom()
        assert_is_instance(atom, (Atom, BasisAtom))

    def test12(self):
        lattice = Crystal2DLattice.square(a=1.0)
        elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
        for element in elements:
            atom = BasisAtom(element, lattice=lattice)
            assert_equal(atom.element, element)

        atom = BasisAtom()
        for c in ('x', 'y', 'z'):
            assert_equal(getattr(atom, c), 0.0)

    def test13(self):
        lattice = Crystal2DLattice.square(a=1.0)
        atom = BasisAtom('C', lattice=lattice)
        assert_equal(atom.element, 'C')
        assert_is_instance(atom.lattice, Crystal2DLattice)

    def test14(self):
        lattice = Crystal3DLattice.cubic(a=5.0)
        basis = BasisAtoms(atoms=['C', 'C'])
        assert_true(basis.lattice is None)
        basis.lattice = lattice
        assert_is_instance(basis.lattice, Crystal3DLattice)

    def test15(self):
        atom = StructureAtom(element='C')
        assert_equal(atom.CN, 0)
        atom.CN = 3
        assert_equal(atom.CN, 3)

    def test16(self):
        atoms1 = StructureAtoms()
        for Z in range(100, 0, -1):
            atoms1.append(StructureAtom(Z=Z))
        atoms1.sort(key=lambda a: a.Z)
        atoms2 = StructureAtoms()
        for Z in range(1, 101):
            atoms2.append(StructureAtom(Z=Z))
        assert_equal(atoms1, atoms2)

    def test17(self):
        atoms = self.atoms
        filtered = atoms.filtered_ids(list(range(5, 26)))
        atoms.filter_ids(list(range(5, 26)))
        assert_equal(filtered, atoms)
        assert_true(atoms.Natoms, 20)

    def test18(self):
        atoms = self.atoms
        filtered1 = atoms.filtered(atoms.coordination_numbers > 1)
        filtered2 = atoms.filtered("coordination_numbers > 1")
        atoms.filter(atoms.coordination_numbers > 1)
        assert_equal(filtered1, atoms)
        assert_equal(filtered1, filtered2)

    def test19(self):
        atoms = self.atoms
        atoms.kNN = 6
        atoms.NNrc = 9.0
        atoms.update_attrs()
        atoms = atoms.filtered((atoms.z >= -5) & (atoms.z <= 5))
        for atom in atoms:
            assert_equal(atom.CN, atoms.kNN)

    def test20(self):
        atoms = self.atoms
        atoms.center_centroid()
        atoms.kNN = 6
        atoms.NNrc = 9.0
        new_atoms = atoms.filtered((atoms.z >= -5) & (atoms.z <= 5))
        assert_equal(atoms.kNN, 6)
        assert_equal(atoms.NNrc, 9.0)
        assert_equal(atoms.kNN, new_atoms.kNN)
        assert_equal(atoms.NNrc, new_atoms.NNrc)

    def test21(self):
        atoms = self.atoms
        assert_is_instance(atoms, StructureAtoms)
        atoms.kNN = 10
        atoms.NNrc = 5
        atoms.update_neighbors(cutoffs=[1.5, 2.0, 3.0])
        # print(atoms.neighbors.tolist())
        # print(atoms.neighbor_distances.tolist())
        # print('first_neighbors:\n{}'.format(atoms.first_neighbors))
        # print('second_neighbors:\n{}'.format(atoms.second_neighbors))
        print('len(first_neighbors): {}'.format(
              [len(atom.first_neighbors) for atom in atoms]))
        print('len(second_neighbors): {}'.format(
              [len(atom.second_neighbors) for atom in atoms]))
        print('len(third_neighbors): {}'.format(
              [len(atom.third_neighbors) for atom in atoms]))
        print('len(4th_neighbors): {}'.format(
              [len(atom.get_nth_nearest_neighbors(4)) for atom in atoms]))

    def test22(self):
        atoms = self.atoms
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.set_pbc()
        atoms.update_attrs()
        for atom in atoms:
            assert_true(np.all(atom.bonds.lengths < 1.5))
            # print('atom: {}, bond.lengths: {}'.format(
            #       atom.id, atom.bonds.lengths))

    def test24(self):
        atoms = self.periodic_table
        atoms.kNN = 3
        atoms_cp = atoms.copy()
        assert_equal(atoms.kNN, atoms_cp.kNN)

    def test25(self):
        atoms = self.atoms
        assert_equal(atoms.Natoms, atoms.ids[-1])

    def test26(self):
        atoms = self.atoms
        assert_true(np.allclose(atoms.coords, atoms.atom_tree.data))
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_attrs()
        assert_equal(len(atoms.nearest_neighbors), atoms.Natoms)
        assert_equal(len(atoms.coordination_numbers), atoms.Natoms)

    def test27(self):
        atoms = self.atoms
        assert_equal(atoms.filtered(atoms.coordination_numbers == 1).Natoms,
                     atoms.coordination_counts[1])
        assert_equal(atoms.filtered(atoms.coordination_numbers == 3).Natoms,
                     atoms.coordination_counts[3])

    def test28(self):
        atoms = self.atoms
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_attrs()
        assert_true(np.allclose(atoms.coordination_numbers,
                                atoms.count_neighbors(2.0)))

    def test29(self):
        atoms = self.atoms
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_attrs()
        # print(atoms.bonds.lengths)
        # print(atoms.neighbor_distances)
        bond_lengths = [bond.length for atom in atoms for bond in atom.bonds]
        assert_true(np.allclose(bond_lengths,
                                atoms.neighbor_distances))

    def test30(self):
        atoms = self.atoms
        assert_true(np.allclose(atoms.volume, atoms.bounds.volume))

    def test31(self):
        atoms = self.atoms
        [assert_equal(atoms.get_atom(i),
                      atoms.get_atom_from_vmd_atom_id(i - 1))
         for i in atoms.ids]

    def test32(self):
        atoms = self.atoms
        assert_equal(len(atoms.x), len(atoms.r))
        assert_is_instance(atoms.r, Vectors)
        assert_is_instance(atoms.x, np.ndarray)

    def test33(self):
        atoms = self.atoms
        atoms.compute_POAVs()
        bonds = list(flatten([atom.bonds for atom in atoms]))
        assert_equal(len(bonds), atoms.coordination_numbers.sum())
        # assert_equal(bonds.Nbonds, atoms.coordination_numbers.sum())
        for atom in atoms:
            if atom.bonds.Nbonds > 1:
                for j, bond in enumerate(atom.bonds):
                    assert_true(np.allclose(bond.vector,
                                            atom.bonds.vectors[j]))
                    assert_equal(bond.length, atom.bonds.lengths[j])

    def test34(self):
        fname = resource_filename('sknano', 'data/nanotubes/0500_5cells.data')
        data_atoms = DATAReader(fname).atoms
        data_atoms.compute_POAVs()
        atoms = self.atoms
        assert_equal(atoms.Natoms, data_atoms.Natoms)

    def test35(self):
        atoms = self.atoms
        atoms.compute_POAVs()
        atoms.filter(atoms.coordination_numbers == 3)
        atom = atoms[10]
        assert_equal(atom.bonds.Nbonds, 3)
        for POAV in ('POAV1', 'POAV2', 'POAVR'):
            atom_POAV = getattr(atom, POAV)

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

    def test36(self):
        atoms = self.atoms
        # atoms.NNrc = 2.0
        atoms.set_pbc()
        atoms.compute_POAVs()

        for atom in atoms:
            # print('atom{}: {}'.format(atom.id, atom))
            for POAV in ('POAV1', 'POAV2', 'POAVR'):
                if getattr(atom, POAV) is not None:
                    atom_POAV = getattr(atom, POAV)
                    sigma_pi_angles = np.degrees(atom_POAV.sigma_pi_angles)
                    assert_false(np.all(np.isclose(sigma_pi_angles,
                                                   3 * [np.nan],
                                                   equal_nan=True)))
                    # print('atom{}.{}.sigma_pi_angles:\n{}'.format(
                    #     atom.id, POAV, sigma_pi_angles))

                    pyramidalization_angles = \
                        np.degrees(atom_POAV.pyramidalization_angles)
                    # print('atom{}.{}.pyramidalization_angles:\n{}'.format(
                    #     atom.id, POAV, pyramidalization_angles))
                    assert_false(np.all(np.isclose(pyramidalization_angles,
                                                   3 * [np.nan],
                                                   equal_nan=True)))

                    misalignment_angles = \
                        np.degrees(atom_POAV.misalignment_angles)
                    # print('atom{}.{}.misalignment_angles:\n{}\n'.format(
                    #     atom.id, POAV, misalignment_angles))
                    assert_false(np.all(np.isclose(misalignment_angles,
                                                   3 * [np.nan],
                                                   equal_nan=True)))

    def test37(self):
        atoms = self.atoms
        atoms.sort(key=attrgetter('CN'))

    def test38(self):
        atoms = self.atoms
        atoms.rotate(axis=[1, 1, 1], angle=np.pi/3)
        # print('inertia_tensor: {}'.format(atoms.inertia_tensor))
        # print('principal_axes: {}'.format(atoms.principal_axes))

    def test39(self):
        atoms = self.atoms
        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, maxlength=20,
                                  retcodes=('ring_counter', 'Rn_counter'))
        assert_equal(ring_cntr[6], 40)
        assert_equal(ring_cntr[10], 9)
        assert_equal(Rn_cntr[0], 10)
        assert_equal(Rn_cntr[2], 10)
        assert_equal(Rn_cntr[3], 10)
        assert_equal(Rn_cntr[4], 70)

    def test40(self):
        atoms = self.atoms
        atoms.set_pbc()
        atoms.update_attrs()
        # print(atoms.nn_adjacency_matrix)
        # print(atoms.nn_adjacency_map)
        # print(atoms.nn_adjacency_list)
        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, maxlength=20,
                                  retcodes=('ring_counter', 'Rn_counter'))
        # print(rings)
        assert_equal(ring_cntr[6], 50)
        assert_equal(ring_cntr[10], 10)
        assert_true(len(list(ring_cntr)) == 2)
        assert_true(set(Rn_cntr), {4})

    def test41(self):
        atoms = self.buckyball
        bonds = atoms.bonds
        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, maxlength=20,
                                  retcodes=('ring_counter', 'Rn_counter'))
        bonds = atoms.bonds
        assert_equal(ring_cntr[6], 20)
        assert_equal(ring_cntr[5], 12)
        assert_equal(ring_cntr[18], 10)
        assert_equal(Rn_cntr[6], 60)
        F = ring_cntr[6] + ring_cntr[5]
        E = bonds.unique_set.Nbonds
        V = atoms.Natoms
        assert_true(V - E + F == 2)

    def test42(self):
        atoms = self.atoms
        atoms.set_pbc()
        atoms.update_attrs()
        atoms.update_neighbor_lists()
        assert_equal(set(atoms.coordination_numbers), {3})
        assert_equal(list(zip(atoms.idx, atoms.nn_idx)), atoms.bonds.indices)

    # def test43(self):
    #     atoms = self.atoms
    #     atoms.set_pbc()
    #     atoms.update_attrs()
    #     G = atoms.graph
    #     print('Graph: {}'.format(list(G.edges())))
    #     print('bfs_edges: {}'.format(list(nx.bfs_edges(G, 0))))
    #     print('len(bfs_edges): {}'.format(len(list(nx.bfs_edges(G, 0)))))
    #     print('bfs_tree: {}'.format(list(nx.bfs_tree(G, 0))))
    #     # print('nx.clustering: {}'.format(nx.clustering(G)))
    #     # print(atoms.nn_adjacency_matrix)
    #     print('bfs_predecessors: {}'.format(dict(nx.bfs_predecessors(G, 0))))
    #     print('bfs_successors: {}'.format(dict(nx.bfs_successors(G, 0))))

    def test43(self):
        atoms = self.atoms
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_attrs()
        for atom in atoms:
            assert_true(np.all(atom.bonds.lengths < 1.5))

if __name__ == '__main__':
    nose.runmodule()
