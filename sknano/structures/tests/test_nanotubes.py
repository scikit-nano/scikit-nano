# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import unittest

from sknano.structures import SWNT, MWNT, SWNTBundle, MWNTBundle


class TestNanotube(unittest.TestCase):

    def setUp(self):
        self.n = 10
        self.m = 10
        self.Ch = (self.n, self.m)

    def test_swnt(self):
        swnt = SWNT(n=self.n, m=self.m)
        self.assertEqual(self.swnt.n, self.n)

        mwnt = MWNT(n=self.n, m=self.m)

if __name__ == '__main__':
    unittest.main()
