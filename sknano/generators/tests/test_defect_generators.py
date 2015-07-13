# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import *

import numpy as np

# from sknano.generators import DefectGenerator
from sknano.io import XYZReader
from sknano.testing import GeneratorTestFixtures


# class TestCase(GeneratorTestFixtures):

#     def test1(self):
#         generator = DefectGenerator()
#         assert_true(isinstance(generator, DefectGenerator))

if __name__ == '__main__':
    nose.runmodule()
