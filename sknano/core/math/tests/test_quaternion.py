#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import warnings

import nose
from nose.tools import *

import numpy as np

from sknano.core.math import Quaternion


def test1():
    q = Quaternion()


if __name__ == '__main__':
    nose.runmodule()
