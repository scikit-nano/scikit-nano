# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename
import unittest

import nose
from nose.tools import *
from sknano.io import DUMPData, DUMPReader, DUMPWriter


#def test_reader():
#    dumpfile = resource_filename('sknano', 'data/lammpstrj/dump.peptide')
#    dumpdata = DUMPReader(dumpfile)

class DUMPTestFixture(unittest.TestCase):

    def setUp(self):
        # dumpfile = \
        #     resource_filename('sknano', 'data/lammpstrj/1010+ion.dump')
        dumpfile = \
            resource_filename('sknano', 'data/lammpstrj/0500_29cells.dump')
        self.dump = \
            DUMPReader(dumpfile, attrmap={'c_peratom_pe': 'pe',
                                          'c_peratom_ke': 'ke'})


class TestCase(DUMPTestFixture):
    def test1(self):
        print(self.dump.dumpattrs)
        print(self.dump.dumpattrs2str())
        print(self.dump.atomattrs)


if __name__ == '__main__':
    nose.runmodule()
