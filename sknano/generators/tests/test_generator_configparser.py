#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import os

import nose
from nose.tools import assert_equal

# import numpy as np

from sknano.testing import GeneratorTestFixture, generate_structure
from sknano.generators import GeneratorConfigParser


class Tests(GeneratorTestFixture):

    def test1(self):
        config = os.path.join(os.path.dirname(__file__),
                              'test_generator_configparser.cfg')
        # print(config)
        parser = GeneratorConfigParser(cfgfile=config)
        # print(parser)

    def test2(self):
        structure = generate_structure(generator_class='SWNTGenerator',
                                       n=10, m=5, nz=5)
        # print(structure)
        parser = GeneratorConfigParser(structure=structure)


if __name__ == '__main__':
    nose.runmodule()
