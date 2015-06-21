# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.structures import MWNTBundle


def test1():
    bundle = MWNTBundle(Ch=[(5, 5), (10, 10)], nx=2, ny=2)
    assert_equal(bundle.Nwalls, 2)
    assert_equal(bundle.Ntubes, 4)


def test2():
    bundle = MWNTBundle(max_walls=5, nx=5, ny=2)
    assert_equal(bundle.Nwalls, 5)
    assert_equal(bundle.Ntubes, 10)


if __name__ == '__main__':
    nose.runmodule()
