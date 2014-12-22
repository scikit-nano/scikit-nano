# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import MWNTBundleGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        bundle = MWNTBundleGenerator(max_shells=2, Lz=1.0,
                                     bundle_geometry='hexagon')
        bundle.save_data()
        self.tmpdata.append(bundle.fname)
        bundle.save_data(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test2(self):
        Ch = [(5,5), (10,10), (15, 15), (20, 20)]
        bundle = MWNTBundleGenerator(Ch=Ch, Lz=0.5,
                                     bundle_geometry='hexagon')
        bundle.save_data()
        self.tmpdata.append(bundle.fname)
        bundle.save_data(structure_format='data')
        self.tmpdata.append(bundle.fname)

    def test3(self):
        Ch = [(3,3), (4, 4), (5,5)]
        bundle = MWNTBundleGenerator(Ch=Ch, nx=5, ny=2, Lz=0.5)
        bundle.save_data()
        self.tmpdata.append(bundle.fname)
        bundle.save_data(structure_format='data')
        self.tmpdata.append(bundle.fname)


if __name__ == '__main__':
    nose.runmodule()
