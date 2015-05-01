# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core import timethis, typeassert, typed_property, \
    ClassSignature


@timethis
def test1():
    class Stock:
        name = typed_property('name', str)
        shares = typed_property('shares', int)
        price = typed_property('price', float)

        def __init__(self, name, shares, price):
            self.name = name
            self.shares = shares
            self.price = price

    s = Stock('ACME', 50, 91.1)
    print(s.name)
    print(s.shares)
    print(s.price)
    assert_true(s.name == 'ACME')
    assert_true(s.shares == 50)
    assert_true(s.price == 91.1)

    with assert_raises(TypeError):
        s.name = 123456

    with assert_raises(TypeError):
        s.shares = '10'

    with assert_raises(TypeError):
        s.price = '-91.1'


def test2():
    @typeassert(int, int)
    def add(x, y):
        return x + y

    assert_equal(add(2, 2), 2 + 2)

    with assert_raises(TypeError):
        add(2, 'hello')


def test3():
    class Stock(ClassSignature):
        _fields = ['name', 'shares', 'price']

    s = Stock('ACME', 50, 91.1)
    print(s.name)
    print(s.shares)
    print(s.price)
    assert_true(s.name == 'ACME')
    assert_true(s.shares == 50)
    assert_true(s.price == 91.1)


if __name__ == '__main__':
    nose.runmodule()
