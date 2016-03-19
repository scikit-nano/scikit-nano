#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal
# import numpy as np

from sknano.core.math import abs_cap


def test1():
    assert_equal(abs_cap(1.1), 1.0)

def test2():
    assert_equal(abs_cap(-1.1), -1.0)


if __name__ == '__main__':
    nose.runmodule()
