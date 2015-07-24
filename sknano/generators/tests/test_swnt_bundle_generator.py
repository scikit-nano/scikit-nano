# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import SWNTGenerator, SWNTBundleGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        bundle = SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
                                     bundle_geometry='hexagon')
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test2(self):
        bundle = SWNTBundleGenerator(n=10, m=0, nx=10, ny=3, nz=5)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test3(self):
        bundle = \
            SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
                                bundle_packing='ccp')
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))
        bundle.save()
        self.tmpdata.append(bundle.fname)
        bundle.save(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test4(self):
        swnt = SWNTGenerator(n=10, m=5, nz=1)
        bundle = \
            SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
                                bundle_packing='ccp')
        assert_true(bundle.Natoms, bundle.Ntubes * swnt.Natoms)
        assert_equal(bundle.Ntubes, len(set(bundle.mol_ids)))


if __name__ == '__main__':
    nose.runmodule()
