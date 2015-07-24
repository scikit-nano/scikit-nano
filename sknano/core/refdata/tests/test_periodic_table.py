# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_names, element_symbols


def test1():
    for data in (atomic_numbers, atomic_number_symbol_map,
                 element_names, element_symbols):
        assert_equal(len(data), len(element_symbols))


def test2():
    print(atomic_masses)
    print(atomic_mass_symbol_map)
    print(list(atomic_masses.keys()))
    assert_equal(element_symbols, list(atomic_masses.keys()))


def test3():
    print(atomic_numbers)
    print(atomic_number_symbol_map)
    print(list(atomic_numbers.keys()))
    assert_equal(element_symbols, list(atomic_numbers.keys()))


def test4():
    assert_equal(len(element_names), len(element_symbols))
    assert_equal(list(atomic_masses.keys()), list(atomic_numbers.keys()))


if __name__ == '__main__':
    nose.runmodule()
