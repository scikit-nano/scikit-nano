# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_almost_equal  # , assert_true
from sknano.structures import SWNT
from sknano.testing import generate_structure


def test1():
    swnt = SWNT(n=10, m=10)
    print(swnt)
    assert_equal(swnt.n, 10)
    assert_equal(swnt.m, 10)
    assert_equal(swnt.element1, 'C')
    assert_equal(swnt.element2, 'C')
    assert_equal(swnt.nz, 1.0)
    assert_equal(swnt.t1, 1)
    assert_equal(swnt.t2, -1)
    assert_equal(swnt.d, 10)
    assert_equal(swnt.dR, 30)
    assert_equal(swnt.N, 20)
    assert_equal(swnt.R, (1, 0))
    assert_almost_equal(swnt.chiral_angle, 30.0)
    assert_almost_equal(swnt.Ch, 42.6, places=2)
    assert_almost_equal(swnt.T, 2.46, places=2)
    assert_almost_equal(swnt.dt, 13.56, places=2)
    assert_almost_equal(swnt.rt, 6.78, places=2)
    assert_equal(swnt.Ntubes, 1)


def test2():
    swnt = SWNT(n=20, m=10)
    print(swnt)
    assert_equal(swnt.element1, 'C')
    assert_equal(swnt.element2, 'C')
    assert_equal(swnt.n, 20)
    assert_equal(swnt.m, 10)
    assert_equal(swnt.nz, 1.0)
    assert_equal(swnt.t1, 4)
    assert_equal(swnt.t2, -5)
    assert_equal(swnt.d, 10)
    assert_equal(swnt.dR, 10)
    assert_equal(swnt.N, 140)
    assert_equal(swnt.R, (1, -1))
    assert_almost_equal(swnt.chiral_angle, 19.11, places=2)
    assert_almost_equal(swnt.Ch, 65.07, places=2)
    assert_almost_equal(swnt.T, 11.27, places=2)
    assert_almost_equal(swnt.dt, 20.71, places=2)
    assert_almost_equal(swnt.rt, 10.36, places=2)
    assert_equal(swnt.Ntubes, 1)


def test3():
    swnt = SWNT(n=20, m=0)
    print(swnt)
    assert_equal(swnt.element1, 'C')
    assert_equal(swnt.element2, 'C')
    assert_equal(swnt.n, 20)
    assert_equal(swnt.m, 0)
    assert_equal(swnt.nz, 1.0)
    assert_equal(swnt.t1, 1)
    assert_equal(swnt.t2, -2)
    assert_equal(swnt.d, 20)
    assert_equal(swnt.dR, 20)
    assert_equal(swnt.N, 40)
    assert_equal(swnt.R, (1, -1))
    assert_almost_equal(swnt.chiral_angle, 0.0, places=2)
    assert_almost_equal(swnt.Ch, 49.2, places=1)
    assert_almost_equal(swnt.T, 4.26, places=2)
    assert_almost_equal(swnt.dt, 15.7, places=1)
    assert_almost_equal(swnt.rt, 7.8, places=1)
    assert_equal(swnt.Ntubes, 1)


def test4():
    swnt = SWNT((10, 5))

    assert_equal(swnt.element1, 'C')
    assert_equal(swnt.element2, 'C')
    print(swnt)
    print(swnt.basis)
    print(swnt.unit_cell)
    print(type(swnt.basis))

    print(swnt)
    print(swnt.basis)
    print(swnt.basis[:2])
    print(type(swnt.basis[:2]))
    print(swnt.crystal_cell)
    print(swnt.crystal_cell.basis)

    swnt.element1 = 'N'
    assert_equal(swnt.element1, 'N')
    print(swnt)

    swnt.element2 = 'Ar'
    assert_equal(swnt.element2, 'Ar')
    print(swnt)
    assert_equal(swnt.unit_cell.basis.symbols.tolist()[:2], ['N', 'Ar'])
    assert_equal(swnt.crystal_cell.basis.symbols.tolist()[:2], ['N', 'Ar'])
    assert_equal(swnt.basis, ['N', 'Ar'])


def test5():
    swnt = SWNT((10, 5), basis=['B', 'N'])
    assert_equal(swnt.element1, 'B')
    assert_equal(swnt.element2, 'N')
    print(swnt.unit_cell.basis.symbols)
    assert_equal(swnt.unit_cell.basis.symbols.tolist()[:2], ['B', 'N'])
    print(swnt)
    print(swnt.unit_cell.basis)

    swnt.element1 = 'N'
    swnt.element2 = 'B'

    assert_equal(swnt.element1, 'N')
    assert_equal(swnt.element2, 'B')

    print(swnt)
    print(swnt.unit_cell)
    assert_equal(swnt.unit_cell.basis.symbols.tolist()[:2], ['N', 'B'])
    assert_equal(swnt.crystal_cell.basis.symbols.tolist()[:2], ['N', 'B'])
    assert_equal(swnt.basis, ['N', 'B'])


def test6():
    swnt = SWNT((10, 5), element1='B', element2='N')
    assert_equal(swnt.element1, 'B')
    assert_equal(swnt.element2, 'N')
    print(swnt)
    print(swnt.unit_cell)
    swnt.element1 = 'N'
    swnt.element2 = 'B'

    assert_equal(swnt.element1, 'N')
    assert_equal(swnt.element2, 'B')

    print(swnt)
    print(swnt.unit_cell)

    assert_equal(swnt.unit_cell.basis.symbols.tolist()[:2], ['N', 'B'])
    assert_equal(swnt.crystal_cell.basis.symbols.tolist()[:2], ['N', 'B'])
    assert_equal(swnt.basis, ['N', 'B'])


def test7():
    structure = generate_structure(generator_class='SWNTGenerator',
                                   n=5, m=0, nz=2)
    print(structure)
    print(structure.crystal_cell)
    print(structure.crystal_cell.lattice)
    print(structure.crystal_cell.basis)
    print(structure.unit_cell)
    print(structure.unit_cell.basis)
    assert_equal(2 * structure.unit_cell.basis.Natoms,
                 structure.crystal_cell.basis.Natoms)

if __name__ == '__main__':
    nose.runmodule()
