# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.structures import Fullerene


def test1():
    fullerene = Fullerene(20)
    assert_equal(fullerene.N, 20)


if __name__ == '__main__':
    nose.runmodule()
