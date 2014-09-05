# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.structures import Graphene


def test1():
    graphene = Graphene(length=5, width=5, edge='AC')
    assert_equals(graphene.length, 5)
    assert_equals(graphene.width, 5)
    assert_equals(graphene.edge, 'AC')


if __name__ == '__main__':
    nose.runmodule()
