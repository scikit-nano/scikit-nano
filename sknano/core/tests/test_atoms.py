#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
import unittest

from .._atom import Atom


class TestAtom(unittest.TestCase):

    def setUp(self):
        self.atom = Atom()

if __name__ == '__main__':
    unittest.main()
