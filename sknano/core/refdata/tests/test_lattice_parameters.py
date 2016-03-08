# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core.refdata import lattice_parameters


def test1():
    for k, v in lattice_parameters.items():
        print('k: {}'.format(k))
        print('v: {}'.format(v))
        if isinstance(v, dict):
            for kk, vv in v.items():
                print('kk: {}'.format(kk))
                print('vv: {}'.format(vv))


if __name__ == '__main__':
    nose.runmodule()
