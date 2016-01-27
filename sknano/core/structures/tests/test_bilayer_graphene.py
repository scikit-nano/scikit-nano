# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

import numpy as np

from sknano.structures import BilayerGraphene


def test1():
    blg = BilayerGraphene(armchair_edge_length=5, zigzag_edge_length=5,
                          layer_rotation_increment=45, degrees=True)
    assert_equal(blg.nlayers, 2)
    assert_equal(blg.layer_rotation_angles[-1], np.pi/4)


if __name__ == '__main__':
    nose.runmodule()
