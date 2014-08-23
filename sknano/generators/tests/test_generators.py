# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from pprint import pprint

import unittest

from sknano.generators import SWNTGenerator, MWNTGenerator, \
    SWNTBundleGenerator, MWNTBundleGenerator


class TestSWNTGenerators(unittest.TestCase):

    def test_swnt_generator(self):
        print('\nTesting SWNTGenerator\n')
        print('SWNTGenerator.__bases__:')
        pprint(SWNTGenerator.__bases__)
        print('\nSWNTGenerator.__mro__:')
        pprint(SWNTGenerator.__mro__)
        swnt = SWNTGenerator(n=10, m=10)
        pprint('test instance: {}'.format(swnt))
        swnt.save_data()
        swnt.save_data(structure_format='data')
        #SWNTGenerator(n=10, m=10).save_data(structure_format='data')

    def test_swntbundle_generator(self):
        print('\nTesting SWNTBundleGenerator\n')
        print('SWNTBundleGenerator.__bases__:')
        pprint(SWNTBundleGenerator.__bases__)
        print('\nSWNTBundleGenerator.__mro__:')
        pprint(SWNTBundleGenerator.__mro__)
        #SWNTBundleGenerator(n=10, m=0, nx=10, ny=3, nz=5).save_data()
        #SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
        #                    bundle_packing='ccp').save_data()
        swntbundle = SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
                                         bundle_geometry='hexagon')
        swntbundle.save_data()
        swntbundle.save_data(structure_format='data')


class TestMWNTGenerators(unittest.TestCase):

    def test_mwnt_generator(self):
        print('\nTesting MWNTGenerator\n')
        print('MWNTGenerator.__bases__:')
        pprint(MWNTGenerator.__bases__)
        print('\nMWNTGenerator.__mro__:')
        pprint(MWNTGenerator.__mro__)
        mwnt = MWNTGenerator(n=20, m=20, max_shells=3, Lz=1.0, fix_Lz=True)
        pprint('test instance: {}'.format(mwnt))
        mwnt.save_data()
        mwnt.save_data(structure_format='data')

    def test_mwntbundle_generator(self):
        print('\nTesting MWNTBundleGenerator\n')
        print('MWNTBundleGenerator.__bases__:')
        pprint(MWNTBundleGenerator.__bases__)
        print('\nMWNTBundleGenerator.__mro__:')
        pprint(MWNTBundleGenerator.__mro__)
        mwntbundle = MWNTBundleGenerator(n=40, m=40, max_shells=5,
                                         Lz=1.0, fix_Lz=True,
                                         bundle_geometry='hexagon')
        mwntbundle.save_data()
        mwntbundle.save_data(structure_format='data')


if __name__ == '__main__':
    unittest.main()
