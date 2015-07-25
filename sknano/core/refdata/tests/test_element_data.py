# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core.refdata import element_data, element_symbols


def test1():
    print(element_data)
    assert_equal(len(element_data.keys()), len(element_symbols))


if __name__ == '__main__':
    nose.runmodule()
