# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core import cyclic_pairs


def test1():
    lst = ['a', 'b', 'c']
    assert_equal(cyclic_pairs(lst), [('a', 'b'), ('b', 'c'), ('c', 'a')])


if __name__ == '__main__':
    nose.runmodule()
