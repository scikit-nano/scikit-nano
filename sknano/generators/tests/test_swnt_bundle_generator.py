# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import SWNTGenerator, MWNTGenerator, \
    SWNTBundleGenerator, MWNTBundleGenerator


def test_swnt_generator():
    swnt = SWNTGenerator(n=10, m=10)
    swnt.save_data()
    swnt.save_data(structure_format='data')


def test_swntbundle_generator():
    #SWNTBundleGenerator(n=10, m=0, nx=10, ny=3, nz=5).save_data()
    #SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
    #                    bundle_packing='ccp').save_data()
    swntbundle = SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
                                     bundle_geometry='hexagon')
    swntbundle.save_data()
    swntbundle.save_data(structure_format='data')


def test_mwnt_generator():
    mwnt = MWNTGenerator(n=20, m=20, max_shells=3, Lz=1.0, fix_Lz=True)
    mwnt.save_data()
    mwnt.save_data(structure_format='data')


def test_mwntbundle_generator():
    mwntbundle = MWNTBundleGenerator(n=40, m=40, max_shells=5, Lz=1.0,
                                     fix_Lz=True, bundle_geometry='hexagon')
    mwntbundle.save_data()
    mwntbundle.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
