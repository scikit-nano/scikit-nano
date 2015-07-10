# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core import ListBasedSet, frozendict


def test1():
    l1 = ListBasedSet('abcdef')
    l2 = ListBasedSet('defghi')
    assert_equal(list(l1 & l2), ['d', 'e', 'f'])


def test2():
    d = {'a': 1, 'b': 2, 'c': 3}
    fd = frozendict(d)
    assert_equal(fd, d)

    d['a'] = fd['c']
    assert_equal(d['a'], fd['c'])

    with assert_raises(TypeError):
        fd['a'] = 1

    with assert_raises(AttributeError):
        fd.update({'a': 1})


if __name__ == '__main__':
    nose.runmodule()
