# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import MWNTBundleGenerator


def test1():
    bundle = MWNTBundleGenerator(max_shells=5, Lz=1.0, fix_Lz=True,
                                 bundle_geometry='hexagon')
    bundle.save_data()
    bundle.save_data(structure_format='data')


def test2():
    Ch = [(10,10), (20,20), (30, 30), (40, 40), (50, 50)]
    bundle = MWNTBundleGenerator(Ch=Ch, Lz=1.0, fix_Lz=True,
                                 bundle_geometry='hexagon')
    bundle.save_data()
    bundle.save_data(structure_format='data')


def test3():
    Ch = [(10,10), (15, 15), (20,20)]
    bundle = MWNTBundleGenerator(Ch=Ch, nx=10, ny=3, Lz=1.0, fix_Lz=True)
    bundle.save_data()
    bundle.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
