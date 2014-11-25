#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename
import numpy as np
from sknano.core.atoms import POAVAtom, POAVAtoms
from sknano.generators import SWNTGenerator
from sknano.io import DATAReader
from sknano.testing import generate_atoms
from six.moves import range
#from sknano.utils.geometric_shapes import Ellipsoid


def test_instantiation():
    from sknano.core.atoms import Atoms, XAtoms, KDTAtoms
    atoms = POAVAtoms()
    assert_is_instance(atoms, (Atoms, XAtoms, KDTAtoms, POAVAtoms))
    swnt_atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=1)
    for atom in swnt_atoms:
        atoms.append(atom)

    assert_equal(len(atoms), len(swnt_atoms))

    atoms = POAVAtoms(atoms=atoms)
    assert_equal(len(atoms), len(swnt_atoms))

    atoms = POAVAtoms(atoms=atoms.data)
    assert_equal(len(atoms), len(swnt_atoms))

    atoms = POAVAtoms(atoms=atoms)
    assert_equal(len(atoms), len(swnt_atoms))


def test_list_methods():
    atoms1 = POAVAtoms()
    for Z in range(100, 0, -1):
        atoms1.append(POAVAtom(Z=Z))
    atoms1.sort(key=lambda a: a.Z)
    atoms2 = POAVAtoms()
    for Z in range(1, 101):
        atoms2.append(POAVAtom(Z=Z))
    assert_equal(atoms1, atoms2)


def test_atom_tree():
    Catoms = POAVAtoms()
    for Z in range(1, 101):
        Catoms.append(POAVAtom(Z=Z))
    assert_equals(Catoms.Natoms, 100)


def test_structure_analysis():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    atoms.assign_unique_ids()
    atoms.compute_POAVs()
    assert_equals(atoms.Natoms, 400)

    atoms = POAVAtoms(atoms=atoms)
    assert_equals(atoms.Natoms, 400)

    atoms.kNN = 3
    atoms.NNrc = 2.0
    NNatoms = atoms.nearest_neighbors
    assert_equals(len(NNatoms), atoms.Natoms)

    atomCN = atoms.coordination_numbers
    assert_equals(len(atomCN), atoms.Natoms)

    a100 = atoms.get_atom(100)
    assert_true(a100 is atoms[99])
    assert_equals(atoms.index(a100), 99)
    a200 = atoms.get_atom(200)
    assert_true(a200 is atoms[199])
    assert_equals(atoms.index(a200), 199)
    a300 = atoms.get_atom(300)
    assert_true(a300 is atoms[299])
    assert_equals(atoms.index(a300), 299)


def test_atom_selections():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=10)
    atoms.assign_unique_ids()
    atoms.compute_POAVs()
    a200 = atoms.get_atom(200)
    assert_true(a200 is atoms[199])
    assert_true(a200 == atoms[199])
    #coords = atoms.coords
    atom_bounds = atoms.bounds
    print('atom bounds: {}'.format(atom_bounds))

    #a200NN = atoms.select_within(Ellipsoid(center=a200.r, r=2.5))
    #assert_equals(a200NN.Natoms, 4)


def test_atom_bonds():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=20, m=10, nz=2)
    atoms.assign_unique_ids()
    atoms.compute_POAVs()
    atoms.NNrc = 2.0
    bonds = atoms.bonds
    #print('bonds: {!s}'.format(bonds))
    assert_equal(len(bonds), atoms.coordination_numbers.sum())
    assert_equal(bonds.Nbonds, atoms.coordination_numbers.sum())

    for i, atom in enumerate(atoms):
        if atom.bonds.Nbonds > 1:
            print('atom.bonds.angles:\n'
                  '{}'.format(np.degrees(atom.bonds.angles)))
            for j, bond in enumerate(atom.bonds):
                assert_true(np.allclose(bond.vector, atom.bonds.vectors[j]))
                assert_equal(bond.length, atom.bonds.lengths[j])


def test_structure_data():
    fname = resource_filename('sknano', 'data/nanotubes/1005_5cells.data')
    swnt = SWNTGenerator(n=10, m=5, nz=5)
    swnt_atoms = swnt.atoms
    swnt_atoms.assign_unique_ids()
    swnt_atoms.compute_POAVs()
    #swnt.save_data(fname=fname, structure_format='data')
    data = DATAReader(fname)
    atoms = data.atoms
    atoms.assign_unique_ids()
    atoms.compute_POAVs()
    assert_equals(swnt_atoms.Natoms, atoms.Natoms)


def test_POAV_angles():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=0, nz=2)
    atoms.assign_unique_ids()
    atoms.update_attrs()
    #atoms.NNrc = 2.0
    #atoms.compute_POAVs()

    for i, atom in enumerate(atoms):
        for POAV in ('POAV1', 'POAV2', 'POAVR'):
            if getattr(atom, POAV) is not None:
                atom_POAV = getattr(atom, POAV)
                print('atom{}.{}.sigma_pi_angles:\n{}'.format(
                    atom.atomID, POAV, np.degrees(atom_POAV.sigma_pi_angles)))
                print('atom{}.{}.pyramidalization_angles:\n{}'.format(
                    atom.atomID, POAV,
                    np.degrees(atom_POAV.pyramidalization_angles)))
                print('atom{}.{}.misalignment_angles:\n{}\n'.format(
                    atom.atomID, POAV,
                    np.degrees(atom_POAV.misalignment_angles)))


#def test_list_mods():
#    atoms = generate_atoms(elements='periodic_table')
#    atoms.assign_unique_ids()
#    #atoms.compute_POAVs()
#    #print(atoms)


if __name__ == '__main__':
    nose.runmodule()
