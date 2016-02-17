#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# from pkg_resources import resource_filename

# from collections import Counter
# from operator import attrgetter

# import networkx as nx
# import numpy as np

import nose
from nose.tools import assert_equal, assert_true

from sknano.testing import AtomsTestFixture


class Tests(AtomsTestFixture):

    def test1(self):
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

    def test2(self):
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

    def test3(self):
        atoms = self.buckyball
        bonds = atoms.bonds

        ring_cntr, Rn_cntr = \
            atoms.analyze_network_py(cutoff=1.5, maxlength=20,
                                     retcodes=('ring_counter', 'Rn_counter'))
        print('Rn_cntr: {}'.format(Rn_cntr))

        ring_cntr, Rn_cntr = \
            atoms.analyze_network(cutoff=1.5, maxlength=20,
                                  retcodes=('ring_counter', 'Rn_counter'))
        print('Rn_cntr: {}'.format(Rn_cntr))
        bonds = atoms.bonds
        assert_equal(ring_cntr[6], 20)
        assert_equal(ring_cntr[5], 12)
        assert_equal(ring_cntr[18], 10)
        assert_equal(Rn_cntr[6], 60)
        F = ring_cntr[6] + ring_cntr[5]
        E = bonds.unique_set.Nbonds
        V = atoms.Natoms
        assert_true(V - E + F == 2)

    def test4(self):
        atoms = self.atoms
        atoms.set_pbc()
        atoms.update_attrs()
        c_rings, c_ring_cntr = \
            atoms.analyze_network(cutoff=1.5, maxlength=50,
                                  retcodes=('rings', 'ring_counter'))
        print(c_ring_cntr)
        assert_true(len(c_rings), sum(c_ring_cntr.values()))

        py_rings, py_ring_cntr = \
            atoms.analyze_network_py(cutoff=1.5, maxlength=50,
                                     retcodes=('rings', 'ring_counter'))
        print(py_ring_cntr)
        assert_true(len(py_rings), sum(py_ring_cntr.values()))

        assert_true(len(c_rings), len(py_rings))
        assert_true(sum(c_ring_cntr.values()), sum(py_ring_cntr.values()))

    # def test5(self):
    #     atoms = self.dumpdata1.trajectory[-1].atoms
    #     # atoms = atoms.filtered((atoms.mol_ids >= 2) & (atoms.mol_ids <= 4))
    #     max_mol_id = 2
    #     atoms = atoms.filtered((atoms.mol_ids <= max_mol_id) &
    #                            (atoms.types == 1))
    #     atoms.verbose = True
    #     assert_equal(atoms.Natoms, 1008 * max_mol_id)
    #     atoms.update_attrs()
    #     rings, ring_cntr = \
    #         atoms.analyze_network(cutoff=2.0, maxlength=40,
    #                               retcodes=('rings', 'ring_counter'))
    #     print(ring_cntr)

    #     rings, ring_cntr = \
    #         atoms.analyze_network_py(cutoff=2.0, maxlength=40,
    #                                  retcodes=('rings', 'ring_counter'))
    #     print(ring_cntr)


if __name__ == '__main__':
    nose.runmodule()
