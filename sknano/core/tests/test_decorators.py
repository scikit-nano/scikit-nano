# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import unittest
#import nose

from sknano.core import deprecated


class TestDecorators(unittest.TestCase):

    @deprecated
    def foo(self):
        return x

    def test_decorator(self):
        self.assertRaises(DeprecationWarning, self.foo(1))


if __name__ == '__main__':
    unittest.main()
