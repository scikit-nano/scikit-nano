# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import SWNTBundleGenerator


def test1():
    bundle = SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
                                 bundle_geometry='hexagon')
    bundle.save_data()
    bundle.save_data(structure_format='data')


def test2():
    bundle = SWNTBundleGenerator(n=10, m=0, nx=10, ny=3, nz=5)
    bundle.save_data()
    bundle.save_data(structure_format='data')


def test3():
    bundle = \
        SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1, bundle_packing='ccp')
    bundle.save_data()
    bundle.save_data(structure_format='data')

if __name__ == '__main__':
    nose.runmodule()
