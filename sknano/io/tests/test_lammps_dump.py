# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal
from sknano.testing import IOTestFixture
# from sknano.io import DUMPReader  # , DUMPData, DUMPWriter


class DUMPTestFixture(IOTestFixture):

    def setUp(self):
        self.dump = self.dump_reader
        self.atoms = self.dump[0].atoms


class Tests(DUMPTestFixture):

    def test1(self):
        atoms = self.atoms
        dump = self.dump
        print('timesteps: {}'.format(dump.timesteps))
        print('Natoms: {}'.format(atoms.Natoms))
        assert_equal(atoms.Natoms, len(atoms.ids))
        print('dump.dumpattrs: {}'.format(dump.dumpattrs))
        print('dump.dumpattrs2str(): {}'.format(dump.dumpattrs2str()))
        print('dump.atomattrs: {}'.format(dump.atomattrs))


if __name__ == '__main__':
    nose.runmodule()
