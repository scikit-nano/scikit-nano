# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from operator import attrgetter

import numpy as np

import nose
from nose.tools import assert_equal, assert_true
from sknano.testing import DUMPTestFixture, GeneratorTestFixture
from sknano.io import DUMPData, DUMPWriter


class Tests(DUMPTestFixture, GeneratorTestFixture):

    def test1(self):
        dump = self.dump_reader
        self.print_dumpattrs(dump)
        atoms = dump[0].atoms
        # print('timesteps: {}'.format(dump.timesteps))
        # print('Natoms: {}'.format(atoms.Natoms))
        assert_equal(atoms.Natoms, len(atoms.ids))

    def test2(self):
        dump = self.dumpdata
        self.print_dumpattrs(dump)
        atoms = dump[0].atoms
        # print('timesteps: {}'.format(dump.timesteps))
        testfile = 'test2.dump'
        self.tmpdata.append(testfile)
        dump.write(testfile)
        test_dump = DUMPData(testfile, dumpattrmap=dump.dumpattrmap,
                             atomattrmap=dump.atomattrmap)
        self.print_dumpattrs(test_dump)
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
                             atomattrmap=dump.atomattrmap)
        self.print_dumpattrs(test_dump)
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
                             atomattrmap=dump.atomattrmap)
        assert_true(np.allclose(dump.timesteps, test_dump.timesteps))
        test_atoms = test_dump[0].atoms
        test_atoms.sort(key=attrgetter('id'))
        # print(test_atoms.ids)
        assert_equal(atoms[:100], test_atoms)

    def test5(self):
        dump = self.dumpdata
        atoms = dump[0].atoms
        atoms.sort(key=attrgetter('id'))
        testfile = 'test5.dump'
        dump.write(testfile, atoms=atoms[:100])
        self.tmpdata.append(testfile)
        test_dump = DUMPData(testfile, dumpattrmap=dump.dumpattrmap,
                             atomattrmap=dump.atomattrmap)
        self.print_dumpattrs(test_dump)
        assert_true(np.allclose(dump.timesteps, test_dump.timesteps))
        test_atoms = test_dump[0].atoms
        test_atoms.sort(key=attrgetter('id'))
        # print(test_atoms.ids)
        assert_equal(atoms[:100], test_atoms)

    def test6(self):
        dump = self.dumpdata
        ss = dump[0]
        atoms = ss.atoms
        atoms.sort(key=attrgetter('id'))
        testfile = 'test6.dump'
        self.tmpdata.append(testfile)
        DUMPWriter.write(testfile, atoms=atoms, timestep=ss.timestep,
                         boxstr=ss.boxstr, domain=ss.domain,
                         style=dump.style,
                         dumpattrs=dump.dumpattrs,
                         dumpattrmap=dump.dumpattrmap,
                         atomattrmap=dump.atomattrmap,
                         center_centroid=False,
                         verbose=True)
        self.print_dumpattrs(dump)
        test_dump = DUMPData(testfile, dumpattrmap=dump.dumpattrmap,
                             atomattrmap=dump.atomattrmap)
        self.print_dumpattrs(test_dump)
        test_atoms = test_dump[0].atoms
        test_atoms.sort(key=attrgetter('id'))
        assert_equal(atoms[:100], test_atoms[:100])

    def test7(self):
        dump = self.dumpdata
        self.print_dumpattrs(dump)
        dump.new_dumpattr('r_mag', values={'r_mag': 'r.magnitude'})
        self.print_dumpattrs(dump)

        ss = dump[0]
        atoms = ss.atoms
        atoms.sort(key=attrgetter('id'))
        testfile = 'test7.dump'
        self.tmpdata.append(testfile)
        DUMPWriter.write(testfile, atoms=atoms, timestep=ss.timestep,
                         boxstr=ss.boxstr, domain=ss.domain,
                         style=dump.style,
                         dumpattrs=dump.dumpattrs,
                         dumpattrmap=dump.dumpattrmap,
                         atomattrmap=dump.atomattrmap,
                         center_centroid=False,
                         verbose=True)
        test_dump = DUMPData(testfile, dumpattrmap=dump.dumpattrmap,
                             atomattrmap=dump.atomattrmap)
        self.print_dumpattrs(test_dump)
        test_atoms = test_dump[0].atoms
        test_atoms.sort(key=attrgetter('id'))
        assert_equal(atoms[:100], test_atoms[:100])

    def test8(self):
        swnt = self.swnt
        # swnt.center_centroid()
        print(swnt.atoms.centroid)
        print(swnt.atoms.lattice.centroid)
        print(swnt.crystal_cell.basis.centroid)
        print(swnt.crystal_cell.lattice.centroid)
        print(swnt.lattice.offset)
        swnt.save(structure_format='dump',
                  dumpattrs=['id', 'type', 'x', 'y', 'z'])
        self.tmpdata.append(swnt.fname)


if __name__ == '__main__':
    nose.runmodule()
