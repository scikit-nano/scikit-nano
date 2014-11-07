# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import GrapheneGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        generator = GrapheneGenerator(width=5, length=5)
        generator.save_data()
        self.tmpdata.append(generator.fname)

    def test2(self):
        generator = GrapheneGenerator(length=20, width=1, edge='ZZ')
        generator.save_data()
        self.tmpdata.append(generator.fname)

    def test3(self):
        generator = GrapheneGenerator(length=20, width=1, edge='AC')
        generator.save_data()
        self.tmpdata.append(generator.fname)


if __name__ == '__main__':
    nose.runmodule()
