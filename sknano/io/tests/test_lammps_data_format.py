# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.io import DATAData, DATAReader, DATAWriter, DATA2XYZConverter


def test_reader():
    infile = '1010r_1cell.data'
    reader = DATAReader(infile)


if __name__ == '__main__':
    nose.runmodule()
