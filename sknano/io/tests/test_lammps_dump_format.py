# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename
import unittest

import nose
from nose.tools import assert_equal
from sknano.io import DUMPReader  # , DUMPData, DUMPWriter


# def test_reader():
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
        self.snapshot0 = self.dump[0]
        self.atoms = self.snapshot0.atoms


class TestCase(DUMPTestFixture):

    def test1(self):
        print('timesteps: {}'.format(self.dump.timesteps))
        print('atoms: {}'.format(self.atoms))
        print('Natoms: {}'.format(self.atoms.Natoms))
        print('atom_ids: {}'.format(self.atoms.ids))
        assert_equal(self.atoms.Natoms, len(self.atoms.ids))
        print(self.dump.dumpattrs)
        print(self.dump.dumpattrs2str())
        print(self.dump.atomattrs)


if __name__ == '__main__':
    nose.runmodule()
