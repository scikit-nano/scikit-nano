# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import SWNTGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        swnt = SWNTGenerator(n=10, m=10)
        swnt.save_data()
        self.tmpdata.append(swnt.fname)
        swnt.save_data(structure_format='data')
        self.tmpdata.append(swnt.fname)

    def test2(self):
        swnt = SWNTGenerator(n=10, m=10, Lz=1.0, fix_Lz=True)
        swnt.save_data()
        self.tmpdata.append(swnt.fname)
        swnt.save_data(structure_format='data')
        self.tmpdata.append(swnt.fname)

if __name__ == '__main__':
    nose.runmodule()
