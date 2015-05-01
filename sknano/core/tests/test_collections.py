# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core import ListBasedSet


def test1():
    l1 = ListBasedSet('abcdef')
    l2 = ListBasedSet('defghi')
    assert_equal(list(l1 & l2), ['d', 'e', 'f'])


if __name__ == '__main__':
    nose.runmodule()
