# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.structures import SWNT, MWNT, SWNTBundle, MWNTBundle


def test_swnt():
    swnt = SWNT(10, 10)
    assertEqual(swnt.n, 10)


def test_swnt_bundle():
    swntbundle = SWNTBundle(10, 10)
    assertEqual(swntbundle.n, 10)


def test_mwnt():
    mwnt = MWNT(10, 10)
    assertEqual(mwnt.n, 10)


def test_mwnt_bundle():
    mwntbundle = MWNTBundle(10, 10)
    assertEqual(mwntbundle.n, 10)

if __name__ == '__main__':
    nose.runmodule()
