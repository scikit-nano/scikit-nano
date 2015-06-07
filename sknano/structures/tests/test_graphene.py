# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.structures import Graphene


def test1():
    graphene = Graphene(armchair_edge_length=5, zigzag_edge_length=5)
    assert_equals(graphene.zigzag_edge_length, 5)
    assert_equals(graphene.armchair_edge_length, 5)
    print(graphene.unit_cell)


if __name__ == '__main__':
    nose.runmodule()
