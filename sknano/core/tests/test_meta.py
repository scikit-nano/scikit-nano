# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core import timethis, Integer, Float, String, \
    UnsignedInteger, UnsignedFloat, SizedString

@timethis
def test1():
    class Stock:
        name = SizedString('name', size=8)
        shares = UnsignedInteger('shares')
        price = UnsignedFloat('price')
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

    with assert_raises(ValueError):
        s.name = '0123456789'

    with assert_raises(TypeError):
        s.name = 123456

    with assert_raises(ValueError):
        s.shares = -10

    with assert_raises(TypeError):
        s.shares = '10'

    with assert_raises(ValueError):
        s.price = -91.1

    with assert_raises(TypeError):
        s.price = '-91.1'

if __name__ == '__main__':
    nose.runmodule()
