# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

from sknano.core.refdata import element_data


def test1():
    print(element_data)


if __name__ == '__main__':
    nose.runmodule()
