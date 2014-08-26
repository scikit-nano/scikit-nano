# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

from sknano.structures import SWNT, MWNT, SWNTBundle, MWNTBundle


class TestSWNT(unittest.TestCase):

    def setUp(self):
        self.swnt = SWNT(10, 10)

    def test_swnt(self):
        self.assertEqual(self.swnt.n, 10)


class TestSWNTBundle(unittest.TestCase):

    def setUp(self):
        self.swntbundle = SWNTBundle(10, 10)

    def test_swnt(self):
        self.assertEqual(self.swntbundle.n, 10)


class TestMWNT(unittest.TestCase):

    def setUp(self):
        self.mwnt = MWNT(10, 10)

    def test_mwnt(self):
        self.assertEqual(self.mwnt.n, 10)


class TestMWNTBundle(unittest.TestCase):

    def setUp(self):
        self.mwntbundle = MWNTBundle(10, 10)

    def test_mwnt(self):
        self.assertEqual(self.mwntbundle.n, 10)

if __name__ == '__main__':
    unittest.main()
