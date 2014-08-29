# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import SWNTGenerator


def test1():
    swnt = SWNTGenerator(n=10, m=10)
    swnt.save_data()
    swnt.save_data(structure_format='data')


def test2():
    swnt = SWNTGenerator(n=10, m=10, Lz=1.0, fix_Lz=True)
    swnt.save_data()
    swnt.save_data(structure_format='data')

if __name__ == '__main__':
    nose.runmodule()
