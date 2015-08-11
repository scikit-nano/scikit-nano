# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_false, assert_true, assert_equal, assert_raises

import logging

from sknano.core import timethis, typeassert, typed_property, \
    ClassSignature, Singleton, Cached, logged
from sknano.core.math import Vector


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


def test4():

    class Spam(metaclass=Singleton):
        pass
        # def __init__(self):
        #     print('Creating spam')

    # print(Spam.__mro__)
    a = Spam()
    b = Spam()
    assert_true(a is b)


def test5():

    class Spam(metaclass=Cached):
        def __init__(self, name):
            print('Creating Spam({!r})'.format(name))
            self.name = name

    print(Spam.__mro__)
    a = Spam('Andrew')
    b = Spam('Python')
    c = Spam('Andrew')
    assert_true(a is c)
    assert_false(a is b)


def test6():
    class NullVector(metaclass=Singleton):
        def __init__(self):
            print('Creating NullVector')
            self.__instance = Vector([0, 0, 0])

        def __repr__(self):
            return repr(self.__instance)

    a = NullVector()
    b = NullVector()
    assert_true(a is b)
    print(a)


def test7():
    @timethis
    @logged
    def add(x, y):
        return x + y

    @logged(level=logging.CRITICAL, name='example')
    @timethis
    def spam():
        print('Spam!')

    logging.basicConfig(level=logging.DEBUG)
    print(add(2, 3))

    add.set_message('Add called')
    print(add(2, 3))

    add.set_level(logging.WARNING)
    print(add(2, 3))


if __name__ == '__main__':
    nose.runmodule()
