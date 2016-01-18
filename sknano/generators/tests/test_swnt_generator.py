# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import assert_equal, assert_true

import numpy as np

from sknano.generators import SWNTGenerator
from sknano.io import XYZReader, DATAReader
from sknano.testing import GeneratorTestFixture


class Tests(GeneratorTestFixture):

    def test1(self):
        swnt = SWNTGenerator(n=5, m=5)
        swnt.save()
        self.tmpdata.append(swnt.fname)
        swnt.save(structure_format='data')
        self.tmpdata.append(swnt.fname)

    def test2(self):
        swnt = SWNTGenerator(n=5, m=5, Lz=1.0, fix_Lz=True)
        swnt.save()
        self.tmpdata.append(swnt.fname)
        swnt.save(structure_format='data')
        self.tmpdata.append(swnt.fname)

    def test3(self):

        new_swnt = SWNTGenerator((10, 10))
        new_swnt.atoms.center_centroid()
        test_swnt = \
            XYZReader(resource_filename('sknano',
                                        'data/nanotubes/1010_1cell.xyz'))
        assert_equal(new_swnt.atoms.Natoms, test_swnt.atoms.Natoms)
        assert_true(np.allclose(new_swnt.atoms.coords, test_swnt.atoms.coords))

    def test4(self):
        bundle = SWNTGenerator(n=5, m=5, nx=3, ny=3, nz=1,
                               bundle_geometry='hexagon')
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test5(self):
        bundle = SWNTGenerator(n=5, m=0, nx=10, ny=3, nz=1)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test6(self):
        bundle = \
            SWNTGenerator(n=5, m=5, nx=3, ny=3, nz=1, bundle_packing='ccp')
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test7(self):
        swnt = SWNTGenerator(n=5, m=5, nz=1)
        bundle = \
            SWNTGenerator(n=5, m=5, nx=3, ny=3, nz=1, bundle_packing='ccp')
        assert_true(bundle.Natoms, bundle.Ntubes * swnt.Natoms)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))

    def test8(self):
        n, m = 5, 5
        nz = 2
        swnt = SWNTGenerator(n=n, m=m, nz=nz)
        assert_equal(swnt.N, 2 * n)
        assert_equal(swnt.Natoms, 2 * swnt.N * swnt.nz)
        assert_equal(swnt.Natoms_per_unit_cell, 2 * swnt.N)
        assert_equal(swnt.nz, nz)
        assert_equal(swnt.Natoms_per_tube, swnt.Natoms_per_unit_cell * swnt.nz)
        atoms = swnt.atoms
        assert_equal(atoms.Natoms, swnt.Natoms_per_tube)

if __name__ == '__main__':
    nose.runmodule()
