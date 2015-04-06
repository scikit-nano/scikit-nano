# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core import cyclic_pairs, dedupe, rezero_array, ListBasedSet


def test1():
    lst = [1, 2, 2, 3, 4, 50, 50, 4, 5]
    assert_equal(list(dedupe(lst)), [1, 2, 3, 4, 50, 5])


def test2():
    lst = ['a', 'b', 'c']
    assert_equal(cyclic_pairs(lst), [('a', 'b'), ('b', 'c'), ('c', 'a')])


def test3():
    l1 = ListBasedSet('abcdef')
    l2 = ListBasedSet('defghi')
    assert_equal(list(l1 & l2), ['d', 'e', 'f'])


if __name__ == '__main__':
    nose.runmodule()
