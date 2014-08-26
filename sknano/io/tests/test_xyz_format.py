# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.io import XYZData, XYZReader, XYZWriter, XYZ2DATAConverter


def test_reader():
    infile = '1010r_1cell.xyz'
    reader = XYZReader(infile)


if __name__ == '__main__':
    nose.runmodule()
