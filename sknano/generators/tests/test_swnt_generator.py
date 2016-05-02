# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import warnings
warnings.simplefilter('always')

import nose
from nose.tools import assert_equal, assert_true

import numpy as np

from sknano.core.structures import SWNT
from sknano.generators import SWNTGenerator
from sknano.io import XYZReader, DATAReader
from sknano.testing import GeneratorTestFixture, generate_structure


class Tests(GeneratorTestFixture):

    def test1(self):
        swnt = SWNTGenerator(n=5, m=5)
        print(swnt)
        swnt.save()
        self.tmpdata.append(swnt.fname)
        swnt.save(structure_format='data')
        self.tmpdata.append(swnt.fname)

    def test2(self):
        swnt = SWNTGenerator(n=5, m=5, L=10.0, fix_L=True)
        print(swnt)
        swnt.save()
        self.tmpdata.append(swnt.fname)
        swnt.save(structure_format='data')
        self.tmpdata.append(swnt.fname)

    def test3(self):

        swnt = SWNTGenerator((10, 10))
        swnt.center_centroid()
        test_swnt = \
            XYZReader(resource_filename('sknano',
                                        'data/nanotubes/1010_1cell.xyz'))
        assert_equal(swnt.Natoms, test_swnt.Natoms)
        assert_true(np.allclose(swnt.coords, test_swnt.coords))

    def test4(self):
        xyz = \
            XYZReader(resource_filename('sknano',
                                        'data/nanotubes/1010_1cell.xyz'))

        data = \
            DATAReader(resource_filename('sknano',
                                         'data/nanotubes/1010_1cell.data'))
        assert_equal(xyz.Natoms, data.Natoms)
        assert_true(np.allclose(xyz.coords, data.coords))

    def test5(self):
        bundle = SWNTGenerator(n=5, m=5, n1=3, n2=3, n3=1,
                               bundle_geometry='hexagon')
        print(bundle)
        print(bundle.scaling_matrix)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        print(bundle.lattice)
        assert_equal(bundle.Ntubes, 7)
        assert_equal(bundle.unit_cell.basis.Natoms, 20)
        assert_equal(bundle.Ntubes * bundle.unit_cell.basis.Natoms,
                     bundle.crystal_cell.basis.Natoms - 40)
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test6(self):
        bundle = SWNTGenerator(n=5, m=0, n1=10, n2=3, n3=1)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test7(self):
        bundle = \
            SWNTGenerator(n=5, m=5, n1=3, n2=3, n3=1, bundle_packing='ccp')
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test8(self):
        swnt = SWNTGenerator(n=5, m=5, n3=1)
        bundle = \
            SWNTGenerator(n=5, m=5, n1=3, n2=3, n3=1, bundle_packing='ccp')
        assert_true(bundle.Natoms, bundle.Ntubes * swnt.Natoms)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))

    def test9(self):
        n, m = 5, 5
        n3 = 2
        swnt = SWNTGenerator(n=n, m=m, n3=n3)
        assert_equal(swnt.N, 2 * n)
        assert_equal(swnt.Natoms, 2 * swnt.N * swnt.n3)
        assert_equal(swnt.Natoms_per_unit_cell, 2 * swnt.N)
        assert_equal(swnt.n3, n3)
        assert_equal(swnt.Natoms_per_tube, swnt.Natoms_per_unit_cell * swnt.n3)
        atoms = swnt.atoms
        assert_equal(atoms.Natoms, swnt.Natoms_per_tube)

    def test10(self):
        structure = generate_structure(generator_class='SWNTGenerator',
                                       n=5, m=0, n3=2)
        assert_equal(2 * structure.unit_cell.basis.Natoms,
                     structure.crystal_cell.basis.Natoms)

    def test11(self):
        swnt = SWNTGenerator((5, 0), L=10, fix_L=True, verbose=False).atoms
        assert_true(swnt.z.max() < 11 and swnt.z.max() > 10)

    def test12(self):
        swnt = SWNT((5, 0), n3=5)
        swnt_generator = SWNTGenerator((5, 0), n3=5).atoms
        assert_equal(swnt.Natoms, swnt_generator.Natoms)

    def test13(self):
        bundle = SWNTGenerator(n=5, m=5, n1=3, n2=3, n3=2)
        # print(bundle.scaling_matrix)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        print(set(bundle.mol_ids))
        print(bundle.lattice)
        assert_equal(bundle.Ntubes, 9)
        assert_equal(bundle.unit_cell.basis.Natoms, 20)
        assert_equal(len(set(bundle.atoms.mol_ids)), bundle.Ntubes)
        assert_equal(bundle.Ntubes * bundle.n3 * bundle.unit_cell.basis.Natoms,
                     bundle.crystal_cell.basis.Natoms)
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test14(self):
        bundle = SWNTGenerator(n=5, m=5, n1=3, n2=3, L=5, fix_L=True)
        print(bundle)
        print(bundle.scaling_matrix)
        # assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        print(bundle.lattice)
        assert_equal(bundle.Ntubes, 9)
        assert_equal(bundle.unit_cell.basis.Natoms, 20)
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

if __name__ == '__main__':
    nose.runmodule()
