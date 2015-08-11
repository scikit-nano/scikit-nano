# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal

from sknano.core import pluralize, ordinal_form


def test1():
    lst = [1, 2, 2, 3, 4, 50, 50, 4, 5]
    assert_equal(pluralize('item', len(lst)), 'items')


def test2():
    lst = [1]
    assert_equal(pluralize('item', len(lst)), 'item')


def test3():
    assert_equal(ordinal_form(0), '0th')
    assert_equal(ordinal_form(1), '1st')
    assert_equal(ordinal_form(2), '2nd')
    assert_equal(ordinal_form(3), '3rd')
    assert_equal(ordinal_form(4), '4th')
    assert_equal(ordinal_form(5), '5th')
    assert_equal(ordinal_form(6), '6th')
    assert_equal(ordinal_form(7), '7th')
    assert_equal(ordinal_form(8), '8th')
    assert_equal(ordinal_form(9), '9th')


def test4():
    for i in range(100):
        print(ordinal_form(i))


if __name__ == '__main__':
    nose.runmodule()
