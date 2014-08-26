# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import MWNTGenerator


def test_mwnt_generator():
    mwnt = MWNTGenerator(n=20, m=20, max_shells=3, Lz=1.0, fix_Lz=True)
    mwnt.save_data()
    mwnt.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
