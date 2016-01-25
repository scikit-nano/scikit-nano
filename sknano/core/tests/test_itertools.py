# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from collections import OrderedDict

import nose
from nose.tools import assert_equal

from sknano.core import cyclic_pairs, take, tabulate, consume, nth, \
    quantify, padnone, ncycles, dotproduct, flatten, repeatfunc, pairwise, \
    grouper, roundrobin, partition, powerset, unique_everseen, \
    unique_justseen, iter_except, first_true, random_product, \
    random_permutation, random_combination, \
    random_combination_with_replacement, ordinal_form


def test1():
    items = ['a', 'b', 'c']
    assert_equal(list(cyclic_pairs(items)),
                 [('a', 'b'), ('b', 'c'), ('c', 'a')])

    items = 'abc'
    assert_equal(list(cyclic_pairs(items)),
                 [('a', 'b'), ('b', 'c'), ('c', 'a')])


def test2():
    items = iter(range(1000))
    assert_equal(take(10, items), list(range(10)))
    assert_equal(take(10, items), list(range(10, 20)))


def test3():
    assert_equal(take(5, tabulate(ordinal_form)),
                 [ordinal_form(i) for i in range(5)])
    assert_equal(take(5, tabulate(ordinal_form, 10)),
                 [ordinal_form(i) for i in range(10, 15)])


def test4():
    i = (x for x in range(10))
    assert_equal(next(i), 0)
    consume(i, 3)
    assert_equal(next(i), 4)


def test5():
    t = tabulate(lambda i: i)
    assert_equal(take(5, pairwise(t)),
                 [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])
    assert_equal(take(5, pairwise(t)),
                 [(6, 7), (7, 8), (8, 9), (9, 10), (10, 11)])
    assert_equal(take(5, pairwise(t)),
                 [(12, 13), (13, 14), (14, 15), (15, 16), (16, 17)])


def test6():
    assert_equal(dotproduct(3 * [1], 3 * [1]), 3)


def test7():
    assert_equal(list(repeatfunc(lambda: 10, 3)), 3 * [10])
    assert_equal(list(repeatfunc(lambda x: x ** 2, 3, 2)), 3 * [4])


def test8():
    assert_equal(list(ncycles("ABC", 3)),
                 ['A', 'B', 'C', 'A', 'B', 'C', 'A', 'B', 'C'])


def test9():
    items = [1, 2, 3, [4, 5, [6, 7], 8], 9]
    assert_equal(list(flatten(items)), list(range(1, 10)))


def test10():
    assert_equal(list(unique_everseen('AAABBBCCDAABBB')),
                 ['A', 'B', 'C', 'D'])
    assert_equal(list(unique_everseen('AAABBCcCAD')),
                 ['A', 'B', 'C', 'c', 'D'])
    assert_equal(list(unique_everseen('AAABBCcCAD', str.lower)),
                 ['A', 'B', 'C', 'D'])


def test11():
    assert_equal(list(unique_justseen('AAABBBCCDAABBB')),
                 ['A', 'B', 'C', 'D', 'A', 'B'])
    assert_equal(list(unique_justseen('AAABBCcCAD', str.lower)),
                 ['A', 'B', 'C', 'A', 'D'])


def test12():
    is_odd = lambda x: x % 2 != 0
    p = partition(is_odd, range(10))
    assert_equal([list(i) for i in p],
                 [list(range(0, 10, 2)), list(range(1, 10, 2))])


def test13():
    assert_equal(list(powerset([1, 2, 3])),
                 [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)])


def test14():
    assert_equal(list(roundrobin('ABC', 'D', 'EF')),
                 ['A', 'D', 'E', 'B', 'F', 'C'])


def test15():
    assert_equal(take(5, padnone(range(3))),
                 [0, 1, 2, None, None])


def test16():
    assert_equal(quantify([True, False, True]), 2)


def test17():
    assert_equal(nth(range(10), 3), 3)
    assert_equal(nth(range(10), 20, 'donkey'), 'donkey')


def test18():
    assert_equal(list(grouper('ABCDEFG', 3, '?')),
                 [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', '?', '?')])


def test19():
    assert_equal(list(iter_except(list(range(3)).pop, IndexError)),
                 [2, 1, 0])

    d = OrderedDict([('a', 5), ('b', 2), ('c', 3)])
    assert_equal(list(iter_except(d.popitem, KeyError)),
                 [('c', 3), ('b', 2), ('a', 5)])


def test20():
    assert_equal(first_true([(), None, '', 0, 1]), 1)
    assert_equal(first_true([(), None, '', 0], 'donkey'), 'donkey')
    is_odd = lambda x: x % 2 != 0
    assert_equal(first_true([2, 4, 6, 8, 11], pred=is_odd), 11)
    assert_equal(first_true([2, 4, 6, 8], default='donkey', pred=is_odd),
                 'donkey')


if __name__ == '__main__':
    nose.runmodule()
