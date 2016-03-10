# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from operator import attrgetter

import numpy as np

import nose
from nose.tools import assert_equal, assert_true
from sknano.testing import IOTestFixture
from sknano.io import DUMPData  # , DUMPWriter


class DUMPTestFixture(IOTestFixture):

    @property
    def dump(self):
        dump = self.dump_reader
        return dump


class Tests(DUMPTestFixture):

    def test1(self):
        dump = self.dump
        print(dump)
        atoms = dump[0].atoms
        # print('timesteps: {}'.format(dump.timesteps))
        # print('Natoms: {}'.format(atoms.Natoms))
        assert_equal(atoms.Natoms, len(atoms.ids))
        print('dump.dumpattrs: {}'.format(dump.dumpattrs))
        print('dump.dumpattrs2str(): {}'.format(dump.dumpattrs2str()))
        print('dump.atomattrs: {}'.format(dump.atomattrs))
        print('dump.atomattrs2str(): {}'.format(dump.atomattrs2str()))

    def test2(self):
        dump = self.dumpdata
        atoms = dump[0].atoms
        # print('timesteps: {}'.format(dump.timesteps))
        testfile = 'test2.dump'
        self.tmpdata.append(testfile)
        dump.write(testfile)
        test_dump = DUMPData(testfile, dumpattrmap=dump.dumpattrmap,
                             elementmap=dump.elementmap)
        assert_true(np.allclose(dump.timesteps, test_dump.timesteps))
        test_atoms = test_dump[0].atoms
        assert_equal(atoms, test_atoms)

    def test3(self):
        dump = self.dumpdata
        ss0 = dump[0]
        atoms = ss0.atoms
        testfile = 'test3.dump.0'
        self.tmpdata.append(testfile)
        dump.write(testfile, snapshot=ss0)
        test_dump = DUMPData(testfile, dumpattrmap=dump.dumpattrmap,
                             elementmap=dump.elementmap)
        assert_equal(dump.timesteps[:1], test_dump.timesteps)
        assert_true(test_dump.Nsnaps == 1)
        test_atoms = test_dump[0].atoms
        assert_equal(atoms, test_atoms)

    def test4(self):
        dump = self.dumpdata
        atoms = dump[0].atoms
        atoms.sort(key=attrgetter('id'))
        testfile = 'test4.dump'
        dump.write(testfile, atoms=atoms[:100])
        self.tmpdata.append(testfile)
        test_dump = DUMPData(testfile, dumpattrmap=dump.dumpattrmap,
                             elementmap=dump.elementmap)
        assert_true(np.allclose(dump.timesteps, test_dump.timesteps))
        test_atoms = test_dump[0].atoms
        test_atoms.sort(key=attrgetter('id'))
        # print(test_atoms.ids)
        assert_equal(atoms[:100], test_atoms)


if __name__ == '__main__':
    nose.runmodule()
