# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import SWNTBundleGenerator


def test1():
    #SWNTBundleGenerator(n=10, m=0, nx=10, ny=3, nz=5).save_data()
    #SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
    #                    bundle_packing='ccp').save_data()
    swntbundle = SWNTBundleGenerator(n=10, m=5, nx=3, ny=3, nz=1,
                                     bundle_geometry='hexagon')
    swntbundle.save_data()
    swntbundle.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
