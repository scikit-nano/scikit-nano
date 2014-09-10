# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import UnrolledSWNTGenerator


def test1():
    unrolled_swnt = UnrolledSWNTGenerator(n=10, m=10)
    unrolled_swnt.save_data()
    unrolled_swnt.save_data(structure_format='data')


def test2():
    unrolled_swnt = UnrolledSWNTGenerator(n=10, m=0)
    unrolled_swnt.save_data()
    unrolled_swnt.save_data(structure_format='data')


def test3():
    unrolled_swnt = UnrolledSWNTGenerator(n=10, m=5)
    unrolled_swnt.save_data()
    unrolled_swnt.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
