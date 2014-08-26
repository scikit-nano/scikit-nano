# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import SWNTGenerator


def test_instance_attributes():
    swnt = SWNTGenerator(n=10, m=10)
    swnt.save_data()
    swnt.save_data(structure_format='data')


if __name__ == '__main__':
    nose.runmodule()
