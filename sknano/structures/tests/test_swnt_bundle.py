# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.structures import SWNTBundle


def test_swnt_bundle():
    swntbundle = SWNTBundle(10, 10)
    assert_equal(swntbundle.n, 10)


if __name__ == '__main__':
    nose.runmodule()
