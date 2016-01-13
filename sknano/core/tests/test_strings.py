# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal

from sknano.core import pluralize, ordinal_form, tuple_expr, list_expr, \
    dict_expr, kwargs_expr, call_signature


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


def test5():
    teststr = "(10, 5)"
    assert_equal(tuple_expr.parseString(teststr)[0], (10.0, 5.0))


def test6():
    teststr = "[5, 5, 5]"
    parsed = list_expr.parseString(teststr, parseAll=True)[0]
    assert_equal(parsed, [5, 5, 5])


def test7():
    teststr = "{'a': 5, 'b': 1, 'c': [1, 1, 2]}"
    assert_equal(dict_expr.parseString(teststr)[0],
                 {'a': 5, 'b': 1, 'c': [1, 1, 2]})


def test8():
    teststr = "dict(a=5, b=1, c=[1, 1, 2])"
    assert_equal(dict_expr.parseString(teststr)[0],
                 {'a': 5, 'b': 1, 'c': [1, 1, 2]})


def test9():
    teststr = "a=5, b=2.5, c=[1, 2, 3]"
    expr = kwargs_expr.parseString(teststr)[0]
    print(expr)
    assert_equal(expr, dict(a=5, b=2.5, c=[1, 2, 3]))


def test10():
    teststr = "(10, 5), 1, 2, 3, a=5, b=10, c=[1, 2, 3]"
    expr = call_signature.parseString(teststr, parseAll=True)[0]
    print(expr)
    args, kwargs = expr
    print('args: {}'.format(args))
    print('kwargs: {}'.format(kwargs))
    assert_equal(args, [(10, 5), 1, 2, 3])
    assert_equal(kwargs, dict(a=5, b=10, c=[1, 2, 3]))

if __name__ == '__main__':
    nose.runmodule()
