# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.generators import FullereneGenerator
from sknano.testing import GeneratorTestFixtures


class TestCase(GeneratorTestFixtures):

    def test1(self):
        buckyball = FullereneGenerator()
        buckyball.save(fname='buckyball.xyz')
        self.tmpdata.append(buckyball.fname)
        buckyball.save(fname='buckyball.data')
        self.tmpdata.append(buckyball.fname)

    def test2(self):
        C70 = FullereneGenerator(N=70)
        C70.save()
        self.tmpdata.append(C70.fname)
        C70.save(structure_format='data')
        self.tmpdata.append(C70.fname)


if __name__ == '__main__':
    nose.runmodule()
