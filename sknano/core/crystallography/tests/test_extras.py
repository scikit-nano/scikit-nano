#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_true
import numpy as np

# from sknano.core import rezero_array
from sknano.core.crystallography import pbc_diff


def test1():
    fc1 = [0.1, 0.1, 0.1]
    fc2 = [0.3, 0.5, 0.9]
    assert_true(np.allclose(pbc_diff(fc1, fc2), [-0.2, -0.4, 0.2]))

    fc3 = [0.9, 0.1, 1.01]
    fc4 = [0.3, 0.5, 0.9]
    assert_true(np.allclose(pbc_diff(fc3, fc4), [-0.4, -0.4, 0.11]))

if __name__ == '__main__':
    nose.runmodule()
