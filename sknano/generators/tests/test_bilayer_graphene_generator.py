# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import BilayerGrapheneGenerator

def test1():
    blg = BilayerGrapheneGenerator(length=1, width=1, layer_rotation_angle=45)



if __name__ == '__main__':
    nose.runmodule()
