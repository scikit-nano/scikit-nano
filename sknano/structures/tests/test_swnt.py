# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.structures import SWNT


def test_swnt():
    swnt = SWNT(10, 10)
    assert_equal(swnt.n, 10)


if __name__ == '__main__':
    nose.runmodule()
