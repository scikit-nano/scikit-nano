# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal
from sknano.core.structures import Fullerene, Fullerenes, load_fullerene_data


def test1():
    for i, (fullerene, data) in enumerate(load_fullerene_data().items()):
        N = Fullerenes.N[i]
        assert_equal(Fullerene(N).N, N)


if __name__ == '__main__':
    nose.runmodule()
