#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

# from collections import Counter
# from operator import attrgetter

# import networkx as nx
import numpy as np

import nose
from nose.tools import assert_equal, assert_true, assert_false, \
    assert_is_instance

from sknano.core import rezero_array
from sknano.core.atoms import Atom, Atoms, BasisAtom, BasisAtoms, \
    StructureAtom, StructureAtoms
from sknano.core.crystallography import Crystal2DLattice, Crystal3DLattice
from sknano.core.geometric_regions import generate_bounding_box
from sknano.core.math import Vector, Vectors, rotation_matrix
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
        assert_equal(atom.q, 0)
        atom.q = 1
        assert_equal(atom.q, 1)

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
        atoms.kNN = 9
        atoms.NNrc = 4.0
        atoms.update_neighbors()
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
        # assert_is_instance(atoms, StructureAtoms)
        atoms.set_pbc('z')
        atoms.kNN = 30
        atoms.NNrc = 1.5
        atoms.update_neighbors(cutoffs=[1.5, 2.0, 3.0])
        assert_true(all([atom.first_neighbors.Natoms == 3 for atom in atoms]))
        assert_true(all([atom.second_neighbors.Natoms == 0 for atom in atoms]))
        assert_true(all([atom.third_neighbors.Natoms == 9 for atom in atoms]))
        # assert_true(all([atom.get_nth_nearest_neighbors(4).Natoms == 6
        #                 for atom in atoms]))

        # print(atoms.neighbors.tolist())
        # print(atoms.neighbor_distances.tolist())
        # print('first_neighbors:\n{}'.format(atoms.first_neighbors))
        # print('second_neighbors:\n{}'.format(atoms.second_neighbors))
        # print('len(first_neighbors): {}'.format(
        #       [len(atom.first_neighbors) for atom in atoms]))
        # print('len(second_neighbors): {}'.format(
        #       [len(atom.second_neighbors) for atom in atoms]))
        # print('len(third_neighbors): {}'.format(
        #       [len(atom.third_neighbors) for atom in atoms]))
        # print('len(4th_neighbors): {}'.format(
        #       [len(atom.get_nth_nearest_neighbors(4)) for atom in atoms]))

    def test22(self):
        atoms = self.atoms
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.set_pbc('xyz')
        atoms.update_neighbors()
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
        atoms.update_neighbors()
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
        atoms.update_neighbors()
        assert_true(np.allclose(atoms.coordination_numbers,
                                atoms.count_neighbors_in_self(2.0)))
        # print(np.sum(atoms.count_neighbors_in_self(2.0)))
        # print(atoms.count_neighbors(atoms.atom_tree, 2.0))
        # print(atoms.query_pairs(2.0))

    def test29(self):
        atoms = self.atoms
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_neighbors()
        # print(atoms.bonds.lengths)
        # print(atoms.neighbor_distances)
        bond_lengths = [bond.length for atom in atoms for bond in atom.bonds]
        assert_true(np.allclose(bond_lengths,
                                atoms.neighbor_distances))

    def test30(self):
        atoms = self.atoms
        lattice = atoms.lattice
        # print(lattice)
        print('lattice.bounding_box: {}'.format(lattice.bounding_box))
        bounding_box = generate_bounding_box(from_lattice=atoms.lattice)
        print(bounding_box)
        assert_equal(bounding_box, lattice.bounding_box)
        assert_true(np.allclose(atoms.bounding_box.volume,
                                bounding_box.volume))

    def test31(self):
        atoms = self.atoms
        [assert_equal(atoms.get_atom(i),
                      atoms.get_atom_from_vmd_index(i - 1))
         for i in atoms.ids]

    def test32(self):
        atoms = self.atoms
        assert_equal(len(atoms.x), len(atoms.r))
        assert_is_instance(atoms.r, Vectors)
        assert_is_instance(atoms.x, np.ndarray)

    def test33(self):
        atoms = self.atoms
        bonds = atoms.all_bonds
        atoms.analyze_POAVs()
        assert_equal(bonds.Nbonds, atoms.coordination_numbers.sum())
        for atom in atoms:
            if atom.bonds.Nbonds > 1:
                for j, bond in enumerate(atom.bonds):
                    assert_true(np.allclose(bond.vector,
                                            atom.bonds.vectors[j]))
                    assert_equal(bond.length, atom.bonds.lengths[j])

    def test34(self):
        fname = resource_filename('sknano', 'data/nanotubes/0500_5cells.data')
        data_atoms = DATAReader(fname).atoms
        data_atoms.analyze_POAVs()
        atoms = self.atoms
        assert_equal(atoms.Natoms, data_atoms.Natoms)

    def test35(self):
        atoms = self.atoms
        atoms.analyze_POAVs()
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
        atoms.set_pbc('xyz')
        atoms.analyze_POAVs()

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
        atoms = self.graphene.atoms
        atoms.set_pbc('xy')
        atoms.kNN = 30
        atoms.NNrc = 1.5
        atoms.update_neighbors(cutoffs=[1.5, 2.0, 2.5, 3.0, 4.0, 4.5])
        assert_true(all([atom.first_neighbors.Natoms == 3 for atom in atoms]))
        assert_true(all([atom.second_neighbors.Natoms == 0 for atom in atoms]))
        assert_true(all([atom.third_neighbors.Natoms == 6 for atom in atoms]))
        assert_true(all([atom.get_nth_nearest_neighbors(4).Natoms == 3
                        for atom in atoms]))
        assert_true(all([atom.get_nth_nearest_neighbors(5).Natoms == 6
                        for atom in atoms]))

    def test38(self):
        atoms = self.atoms
        atoms.rotate(axis=[1, 1, 1], angle=np.pi/3)
        # print('inertia_tensor: {}'.format(atoms.inertia_tensor))
        # print('principal_axes: {}'.format(atoms.principal_axes))

    def test39(self):
        atoms = self.atoms
        atoms.set_pbc('xyz')
        atoms.update_neighbors()
        atoms.update_neighbor_lists()
        assert_equal(set(atoms.coordination_numbers), {3})
        # print(list(zip(atoms.idx, atoms.nn_idx)))
        # assert_equal(list(zip(atoms.idx, atoms.nn_idx)), atoms.bonds.indices)

    def test40(self):
        atoms = self.atoms
        atoms.kNN = 3
        atoms.NNrc = 2.0
        atoms.update_neighbors()
        for atom in atoms:
            assert_true(np.all(atom.bonds.lengths < 1.5))

    def test41(self):
        atoms = self.atoms
        atoms.analyze_POAVs()
        Vpi_vectors, Vpi_vectors_atom_ids = atoms.get_POAV_attr('POAVR', 'Vpi')
        assert_is_instance(Vpi_vectors, Vectors)
        assert_is_instance(Vpi_vectors[0], Vector)
        Vpi_y_components, Vpi_y_components_atom_ids = \
            atoms.get_POAV_attr('POAVR', 'Vpi.y')
        assert_is_instance(Vpi_y_components[0], float)
        assert_equal(len(Vpi_vectors), len(Vpi_y_components))

    def test42(self):
        atoms = self.atoms
        assert_equal(atoms.get_vmd_selection_string('index'),
                     ' '.join(('index', ' '.join(map(str, atoms.indices)))))

    def test43(self):
        atoms = self.dumpdata1.trajectory[-1].atoms
        layer4 = atoms.filtered(atoms.mol_ids == 4)
        assert_equal(layer4.Natoms, (atoms.Natoms - 2) // 10)

    def test44(self):
        atoms = self.graphene.atoms
        print(atoms.Natoms)
        print(atoms.neighbor_cutoffs)
        atoms.set_pbc('xy')
        assert_is_instance(atoms, StructureAtoms)
        atoms.kNN = 12
        # atoms.NNrc = 5
        # atoms.update_neighbors(cutoffs=[1.5, 2.0, 3.0])
        atoms.update_neighbors(cutoffs=[2.0, 3.0])
        atom30 = atoms.get_atom(30)
        assert_equal(atom30.first_neighbors.Natoms, 3)
        print(atom30.first_neighbors.pbc)
        # assert_true(all(atom30.first_neighbors.pbc))
        # print(atom30.neighbor_map)
        print(atom30.second_neighbors.Natoms)
        print(atom30.neighbor_map['2nd'].Natoms)

    def test45(self):
        atoms = self.atoms
        atoms.center_com()
        atoms.rotate(axis='x', angle=np.pi/4)
        R = rotation_matrix(axis='x', angle=-np.pi/4)
        print(R)
        # print(atoms.bounding_box)
        i = atoms.inertia_tensor
        print('inertia_tensor:\n{}'.format(i))
        p = atoms.principal_axes
        print('principal_axes:\n{}'.format(p))
        print('principal_axes.T:\n{}'.format(p.T))
        atoms.rotate(transform_matrix=p.A)
        print(atoms.inertia_tensor)
        # print(p[0])
        # print(p.x)
        w = atoms.principal_moments_of_inertia
        print('principal_moments_of_inertia:\n{}'.format(w))
        W = np.diag(w)
        print('principal_moments_matrix:\n{}'.format(W))
        I = rezero_array(np.asmatrix(p.A) * np.asmatrix(i) * np.asmatrix(p.T),
                         epsilon=1e-10)
        print('rotated inertia_tensor:\n{}'.format(I))

    def test46(self):
        atoms = self.atoms
        bounding_region = atoms.bounding_region
        lattice_region = atoms.lattice_region
        assert_equal(bounding_region, lattice_region)
        bounding_box = atoms.bounding_box
        assert_equal(bounding_box, lattice_region.bounding_box)
        coordinates_bbox = atoms.coordinates_bounding_box
        rminmax = atoms.r.minmax
        print(rminmax)
        assert_equal(coordinates_bbox, coordinates_bbox.bounding_box)


if __name__ == '__main__':
    nose.runmodule()
