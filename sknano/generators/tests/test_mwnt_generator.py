# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import MWNTGenerator


def test1():
    mwnt = MWNTGenerator(max_shells=3, Lz=1.0, fix_Lz=True)
    mwnt.save_data()
    mwnt.save_data(structure_format='data')


def test2():
    mwnt = MWNTGenerator(Ch=[(10, 10), (20, 20), (30, 30), (40, 40), (50, 50)],
                         Lz=1.0, fix_Lz=True)
    mwnt.save_data()
    mwnt.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
