# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import UnrolledSWNTGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        unrolled_swnt = UnrolledSWNTGenerator(n=10, m=10)
        unrolled_swnt.save()
        self.tmpdata.append(unrolled_swnt.fname)
        unrolled_swnt.save(structure_format='data')
        self.tmpdata.append(unrolled_swnt.fname)

    def test2(self):
        unrolled_swnt = UnrolledSWNTGenerator(n=10, m=0)
        unrolled_swnt.save()
        self.tmpdata.append(unrolled_swnt.fname)
        unrolled_swnt.save(structure_format='data')
        self.tmpdata.append(unrolled_swnt.fname)

    def test3(self):
        unrolled_swnt = UnrolledSWNTGenerator(n=10, m=5)
        unrolled_swnt.save()
        self.tmpdata.append(unrolled_swnt.fname)
        unrolled_swnt.save(structure_format='data')
        self.tmpdata.append(unrolled_swnt.fname)


if __name__ == '__main__':
    nose.runmodule()
