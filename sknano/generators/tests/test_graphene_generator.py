# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import GrapheneGenerator


def test1():
    generator = GrapheneGenerator(width=5, length=5)
    generator.save_data()


def test2():
    generator = GrapheneGenerator(length=20, width=1, edge='ZZ')
    generator.save_data()


def test3():
    generator = GrapheneGenerator(length=20, width=1, edge='AC')
    generator.save_data()


if __name__ == '__main__':
    nose.runmodule()
