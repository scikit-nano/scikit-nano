# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core import pluralize


def test1():
    lst = [1, 2, 2, 3, 4, 50, 50, 4, 5]
    assert_equal(pluralize('item', len(lst)), 'items')


def test1():
    lst = [1]
    assert_equal(pluralize('item', len(lst)), 'item')


if __name__ == '__main__':
    nose.runmodule()
