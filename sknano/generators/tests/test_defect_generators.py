# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename

import nose
from nose.tools import assert_raises, assert_true, assert_equal

import numpy as np

from sknano.generators import CrossLinkedDefectGenerator, VacancyGenerator
from sknano.io import XYZReader
from sknano.testing import GeneratorTestFixture, IOTestFixture


class Tests(GeneratorTestFixture, IOTestFixture):

    def test1(self):
        generator = CrossLinkedDefectGenerator()
        assert_true(isinstance(generator, CrossLinkedDefectGenerator))
        # print(generator)

    def test2(self):
        structure = self.swnt.structure
        generator = CrossLinkedDefectGenerator(self.swnt)
        assert_true(isinstance(generator, CrossLinkedDefectGenerator))
        assert_equal(generator.structure, structure)
        # print(generator)

    def test3(self):
        structure = self.xyzdata.structure
        generator = CrossLinkedDefectGenerator(self.xyzdata)
        assert_true(isinstance(generator, CrossLinkedDefectGenerator))
        assert_equal(generator.structure, structure)
        # print(generator)

    def test4(self):
        structure = self.swnt.structure
        with assert_raises(TypeError):
            VacancyGenerator(structure)

    def test5(self):
        structure = self.swnt.structure
        generator = VacancyGenerator(structure, Nsites=5)
        assert_equal(generator.Nsites, 5)

    def test6(self):
        structure = self.swnt.structure
        generator = VacancyGenerator(structure, Nsites=5, Nvacs_per_site=3)
        assert_equal(generator.Nsites, 5)
        assert_equal(generator.Nvacs_per_site, 3)


if __name__ == '__main__':
    nose.runmodule()
