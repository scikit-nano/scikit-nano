# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import MWNTBundleGenerator


def test1():
    bundle = MWNTBundleGenerator(max_shells=2, Lz=1.0, fix_Lz=True,
                                 bundle_geometry='hexagon')
    bundle.save_data()
    bundle.save_data(structure_format='data')


def test2():
    Ch = [(5,5), (10,10), (15, 15), (20, 20)]
    bundle = MWNTBundleGenerator(Ch=Ch, Lz=0.5, fix_Lz=True,
                                 bundle_geometry='hexagon')
    bundle.save_data()
    bundle.save_data(structure_format='data')


def test3():
    Ch = [(3,3), (4, 4), (5,5)]
    bundle = MWNTBundleGenerator(Ch=Ch, nx=5, ny=2, Lz=0.5, fix_Lz=True)
    #bundle.save_data()
    #bundle.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
