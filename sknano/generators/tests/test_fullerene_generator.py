# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import FullereneGenerator


def test1():
    buckyball = FullereneGenerator()
    buckyball.save_data(fname='C60.data')


def test2():
    C70 = FullereneGenerator(N=70)
    C70.save_data()
    C70.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
