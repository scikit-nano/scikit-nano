# -*- coding: utf-8 -*-
"""
========================================================================
Test fixtures (:mod:`sknano.testing._tools`)
========================================================================

.. currentmodule:: sknano.testing._tools

"""
from __future__ import absolute_import, division, print_function

__all__ = ['GeneratorTestFixtures']


import os
import unittest


class GeneratorTestFixtures(unittest.TestCase):

    def setUp(self):
        self.tmpdata = []

    def tearDown(self):
        for f in self.tmpdata:
            try:
                os.remove(f)
            except IOError:
                continue