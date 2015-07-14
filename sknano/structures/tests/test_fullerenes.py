# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.structures import Fullerene, Fullerenes, load_fullerene_data


def test1():
    fullerene = Fullerene(20)
    print(fullerene)
    assert_equal(fullerene.N, 20)


def test2():
    fullerenes = Fullerenes()
    print(fullerenes)


def test3():
    data = load_fullerene_data()
    print(list(data.keys()))
    print(list(data.values()))


if __name__ == '__main__':
    nose.runmodule()
