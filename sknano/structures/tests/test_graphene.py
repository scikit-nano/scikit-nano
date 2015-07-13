# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

import numpy as np

from sknano.structures import Graphene, PrimitiveCellGraphene, \
    ConventionalCellGraphene, GraphenePrimitiveCell, GrapheneConventionalCell


def test1():
    s = Graphene(armchair_edge_length=5, zigzag_edge_length=5)
    assert_equals(s.zigzag_edge_length, 5)
    assert_equals(s.armchair_edge_length, 5)
    assert_true(isinstance(s, ConventionalCellGraphene))
    assert_true(isinstance(s.unit_cell, GrapheneConventionalCell))
    print(s.unit_cell)


def test2():
    s = PrimitiveCellGraphene(edge_length=5)
    assert_equals(s.edge_length, 5)
    assert_true(isinstance(s, PrimitiveCellGraphene))
    assert_true(isinstance(s.unit_cell, GraphenePrimitiveCell))
    print(np.degrees(s.r1.angle(s.r2)))
    print(s.unit_cell)
    print(s.area)
    print(s)


def test3():
    s = ConventionalCellGraphene(armchair_edge_length=5, zigzag_edge_length=5)
    assert_equals(s.zigzag_edge_length, 5)
    assert_equals(s.armchair_edge_length, 5)
    assert_true(isinstance(s, ConventionalCellGraphene))
    assert_true(isinstance(s.unit_cell, GrapheneConventionalCell))
    print(s.unit_cell)
    print(s.area)
    print(s)


def test4():
    s = Graphene.from_conventional_cell(armchair_edge_length=5,
                                        zigzag_edge_length=5)
    assert_equals(s.zigzag_edge_length, 5)
    assert_equals(s.armchair_edge_length, 5)
    assert_true(isinstance(s.unit_cell, GrapheneConventionalCell))
    print(s.unit_cell)
    assert_true(isinstance(s, ConventionalCellGraphene))


def test5():
    s = Graphene.from_primitive_cell(edge_length=5)
    assert_true(isinstance(s, PrimitiveCellGraphene))
    assert_true(isinstance(s.unit_cell, GraphenePrimitiveCell))


if __name__ == '__main__':
    nose.runmodule()
