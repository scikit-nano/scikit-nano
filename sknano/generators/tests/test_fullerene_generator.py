# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from sknano.generators import FullereneGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        buckyball = FullereneGenerator()
        buckyball.save_data(fname='buckyball.xyz')
        self.tmpdata.append(buckyball.fname)
        buckyball.save_data(fname='buckyball.data')
        self.tmpdata.append(buckyball.fname)

    def test2(self):
        C70 = FullereneGenerator(N=70)
        C70.save_data()
        self.tmpdata.append(C70.fname)
        C70.save_data(structure_format='data')
        self.tmpdata.append(C70.fname)


if __name__ == '__main__':
    nose.runmodule()
