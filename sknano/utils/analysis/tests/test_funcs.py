#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import assert_equal

# import numpy as np

from sknano.core import flatten
from sknano.core.atoms import generate_vmd_selection_string
from sknano.core.geometric_regions import Cylinder
from sknano.io import DUMPReader
from sknano.utils.analysis import find_defect_chains, find_target_atom
from sknano.testing import generate_atoms


def test_find_target_atom():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=3)
    atoms.center_centroid()
    # atoms.update_attrs()
    target_atom = find_target_atom(atoms, target_coords=[0.0, 5.0, 0.0],
                                   search_radius=1.0, nearest_target=True)

    print('atoms.z.min(), atoms.z.max(): {}, {}'.format(
        atoms.z.min(), atoms.z.max()))
    print('atom_ids: {}'.format(atoms.atom_ids))
    print('CNs: {}'.format(atoms.coordination_numbers))

    print('target_atom: {}'.format(target_atom))


def test_find_defect_chains():
    dumpfile = \
        resource_filename('sknano',
                          'data/lammpstrj/' +
                          'irradiated_graphene.system.dump.02000')

    dump = DUMPReader(dumpfile, dumpattrmap={'c_atom_pe': 'pe',
                                             'c_atom_ke': 'ke'},
                      atomattrmap={('type', 'element'): {1: 'C', 2: 'Ar'}})
    atoms = dump.trajectory[-1].atoms
    max_mol_id = 3
    atoms = atoms.filtered((atoms.types == 1) & (atoms.mol_ids <= max_mol_id))
    # print(atoms.ids)
    assert_equal(atoms.Natoms, 1008 * max_mol_id)
    cylinder = Cylinder(p1=[0, -25, 0], p2=[0, 25, 0], r=15)
    atoms = atoms.within_region(cylinder)
    atoms.verbose = True
    rings, ring_cntr, Rn_cntr = \
        atoms.analyze_network(cutoff=2.0, max_ring_size=20,
                              retcodes=('rings', 'ring_counter',
                                        'Rn_counter'))

    print('ring_counter: {}'.format(ring_cntr))
    print('Rn_counter: {}'.format(Rn_cntr))
    defect_condition = lambda atom: atom.CN < 3

    all_ring_atoms = []
    # all_ring_defects = []
    for n, list_of_atoms in rings.items():
        if n != 6:
            all_ring_atoms.extend([ring_atom for ring_atoms in list_of_atoms
                                   for ring_atom in ring_atoms])

    for n, list_of_atoms in rings.items():
        if n != 6:
            print("{}-Rings:".format(n))
            print("Ring Counts: {}".format(ring_cntr[n]))
            # print("Rn Counts: {}".format(Rn_cntr[n]))
            atom_indices = [ring_atoms.vmd_indices.tolist() for ring_atoms in
                            list_of_atoms]
            print("Ring IDs:\n{}".format(atom_indices))
            print("VMD selection string:\n{}\n".format(
                  generate_vmd_selection_string(
                      'index', list(set(flatten(atom_indices))))))

            ring_defects = []
            for ring_atoms in list_of_atoms:
                defects = \
                    find_defect_chains(ring_atoms,
                                       defect_condition=defect_condition,
                                       max_length=10,
                                       exclude=all_ring_atoms)
                if len(defects) != 0:
                    ring_defects.extend(defects)

            if len(ring_defects) != 0:
                defect_atom_indices = [defect_atoms.vmd_indices.tolist()
                                       for defect_atoms in ring_defects]
                print("Defect IDs:\n{}".format(defect_atom_indices))
                print("VMD selection string:\n{}\n".format(
                      generate_vmd_selection_string(
                          'index', list(set(flatten(defect_atom_indices))))))

if __name__ == '__main__':
    nose.runmodule()
