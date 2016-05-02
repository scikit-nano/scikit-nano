# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import *

import numpy as np

from sknano.generators import CrossLinkedDefectGenerator
from sknano.io import XYZReader
from sknano.testing import GeneratorTestFixture, IOTestFixture


class Tests(GeneratorTestFixture, IOTestFixture):

    def test1(self):
        generator = CrossLinkedDefectGenerator()
        assert_true(isinstance(generator, CrossLinkedDefectGenerator))
        assert_true(generator.structure is None)
        print(generator)

    def test2(self):
        structure = self.swnt.structure
        generator = CrossLinkedDefectGenerator(self.swnt)
        assert_true(isinstance(generator, CrossLinkedDefectGenerator))
        assert_equal(generator.structure, structure)
        print(generator)

    def test3(self):
        structure = self.xyzdata.structure
        generator = CrossLinkedDefectGenerator(self.xyzdata)
        assert_true(isinstance(generator, CrossLinkedDefectGenerator))
        assert_equal(generator.structure, structure)
        print(generator)


if __name__ == '__main__':
    nose.runmodule()
