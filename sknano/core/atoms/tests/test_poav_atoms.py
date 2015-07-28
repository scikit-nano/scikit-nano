#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from pkg_resources import resource_filename
import numpy as np
from sknano.core.atoms import POAVAtom, POAVAtoms
from sknano.generators import SWNTGenerator
from sknano.io import DATAReader
from sknano.testing import generate_atoms


def test_list_methods():
    atoms1 = POAVAtoms()
    for Z in range(100, 0, -1):
        atoms1.append(POAVAtom(Z=Z))
    atoms1.sort(key=lambda a: a.Z)
    atoms2 = POAVAtoms()
    for Z in range(1, 101):
        atoms2.append(POAVAtom(Z=Z))
    assert_equal(atoms1, atoms2)


def test_atom_bonds():
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
                assert_true(np.allclose(bond.vector, atom.bonds.vectors[j]))
                assert_equal(bond.length, atom.bonds.lengths[j])


def test_structure_data():
    fname = resource_filename('sknano', 'data/nanotubes/1005_5cells.data')
    swnt = SWNTGenerator(n=10, m=5, nz=5)
    swnt_atoms = swnt.atoms
    swnt_atoms.compute_POAVs()
    data = DATAReader(fname)
    atoms = data.atoms
    atoms.compute_POAVs()
    assert_equals(swnt_atoms.Natoms, atoms.Natoms)


def test_POAV_angles():
    atoms = \
        generate_atoms(generator_class='SWNTGenerator', n=10, m=0, nz=2)
    # atoms.NNrc = 2.0
    atoms.compute_POAVs()

    for i, atom in enumerate(atoms):
        print('atom{}: {}'.format(atom.id, atom))
        for POAV in ('POAV1', 'POAV2', 'POAVR'):
            if getattr(atom, POAV) is not None:
                atom_POAV = getattr(atom, POAV)
                print('atom{}.{}.sigma_pi_angles:\n{}'.format(
                    atom.id, POAV, np.degrees(atom_POAV.sigma_pi_angles)))
                print('atom{}.{}.pyramidalization_angles:\n{}'.format(
                    atom.id, POAV,
                    np.degrees(atom_POAV.pyramidalization_angles)))
                print('atom{}.{}.misalignment_angles:\n{}\n'.format(
                    atom.id, POAV,
                    np.degrees(atom_POAV.misalignment_angles)))


if __name__ == '__main__':
    nose.runmodule()
