# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import MWNTBundleGenerator


def test1():
    bundle = MWNTBundleGenerator(n=40, m=40, max_shells=5, Lz=1.0,
                                 fix_Lz=True, bundle_geometry='hexagon')
    bundle.save_data()
    bundle.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
