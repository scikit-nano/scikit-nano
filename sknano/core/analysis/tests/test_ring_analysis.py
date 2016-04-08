#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# from pkg_resources import resource_filename

# from collections import Counter
# from operator import attrgetter

# import numpy as np

import nose
from nose.tools import assert_equal, assert_true

from sknano.core import flatten
from sknano.core.geometric_regions import Cylinder
from sknano.testing import AtomsTestFixture


class Tests(AtomsTestFixture):

    def test1(self):
        atoms = self.atoms
        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, max_ring_size=20,
                                  retcodes=('ring_counter', 'Rn_counter'))
        assert_equal(ring_cntr[6], 40)
        assert_equal(ring_cntr[10], 9)
        assert_equal(Rn_cntr[0], 10)
        assert_equal(Rn_cntr[2], 10)
        assert_equal(Rn_cntr[3], 10)
        assert_equal(Rn_cntr[4], 70)

    def test2(self):
        atoms = self.atoms
        atoms.set_pbc('xyz')
        atoms.update_neighbors()
        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, max_ring_size=20,
                                  retcodes=('ring_counter', 'Rn_counter'))
        # print(rings)
        assert_equal(ring_cntr[6], 50)
        assert_equal(ring_cntr[10], 10)
        assert_true(len(list(ring_cntr)) == 2)
        assert_true(set(Rn_cntr), {4})

    def test3(self):
        atoms = self.buckyball
        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, max_ring_size=20,
                                  retcodes=('ring_counter', 'Rn_counter'))
        bonds = atoms.bonds
        assert_equal(ring_cntr[6], 20)
        assert_equal(ring_cntr[5], 12)
        assert_equal(ring_cntr[18], 10)
        assert_equal(Rn_cntr[6], 60)
        F = ring_cntr[6] + ring_cntr[5]
        E = bonds.unique.Nbonds
        V = atoms.Natoms
        assert_true(V - E + F == 2)

        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, max_ring_size=20,
                                  retcodes=('ring_counter', 'Rn_counter'),
                                  pyversion=True)
        assert_equal(ring_cntr[6], 20)
        assert_equal(ring_cntr[5], 12)
        assert_equal(ring_cntr[18], 10)
        assert_equal(Rn_cntr[6], 60)

    def test4(self):
        atoms = self.atoms
        atoms.set_pbc('xyz')
        atoms.update_neighbors()
        c_rings, c_ring_cntr = \
            atoms.analyze_network(cutoff=1.5, max_ring_size=50,
                                  retcodes=('rings', 'ring_counter'))
        # print(c_ring_cntr)
        assert_true(len(c_rings), sum(c_ring_cntr.values()))
        assert_equal(c_ring_cntr[6], 50)
        assert_equal(c_ring_cntr[10], 10)

        py_rings, py_ring_cntr = \
            atoms.analyze_network(cutoff=1.5, max_ring_size=50,
                                  retcodes=('rings', 'ring_counter'),
                                  pyversion=True)
        # print(py_ring_cntr)
        assert_true(len(py_rings), sum(py_ring_cntr.values()))
        assert_equal(py_ring_cntr[6], 50)
        assert_equal(py_ring_cntr[10], 10)

        assert_true(len(c_rings), len(py_rings))
        assert_true(sum(c_ring_cntr.values()), sum(py_ring_cntr.values()))

    def test5(self):
        reference_angle = 2.094354
        reference_bond = 1.44566

        cylinder = Cylinder(p1=[0, -15, 0], p2=[0, 25, 0], r=10)
        atoms = self.dumpdata1.trajectory[-1].atoms
        atoms = atoms.within_region(cylinder)
        print(atoms.bounding_box)
        atoms.update_neighbors(cutoff=2.0)
        selections = []
        for i in range(1, 9):
            selections.append(atoms.select('molid {}'.format(i)))

        assert_equal(list(range(1, 9)),
                     list(flatten([set(selection.mol_ids)
                                   for selection in selections])))

        # max_mol_id = 2
        # atoms = atoms.filtered((atoms.mol_ids <= max_mol_id) &
        #                        (atoms.types == 1))
        atoms.verbose = True
        # assert_equal(atoms.Natoms, 1008 * max_mol_id)
        atoms.update_neighbors()
        rings, ring_cntr = \
            atoms.analyze_network(cutoff=2.0, max_ring_size=20,
                                  retcodes=('rings', 'ring_counter'))
        print(ring_cntr)
        # atoms.angles_in_degrees = True
        atoms.update_ring_stats(angle_reference=reference_angle,
                                bond_reference=reference_bond)
        # angles = atoms.angles
        # angles.compute_strains(reference_angle)

        # bonds = atoms.bonds
        # bonds.compute_strains(reference_bond)
        print(atoms.ring_stats)

        # rings, ring_cntr = \
        #     atoms.analyze_network(cutoff=2.0, max_ring_size=20,
        #                           retcodes=('rings', 'ring_counter'),
        #                           pyversion=True)
        # print(ring_cntr)


if __name__ == '__main__':
    nose.runmodule()
