#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_true, assert_equal, assert_is_instance
from sknano.core.atoms import BasisAtom, BasisAtoms
from sknano.core.crystallography import Crystal2DLattice, Crystal3DLattice
from sknano.testing import generate_atoms
# from sknano.core.geometric_regions import Ellipsoid


def test1():
    from sknano.core.atoms import Atom
    xatom = BasisAtom()
    assert_is_instance(xatom, (Atom, BasisAtom))


def test2():
    lattice = Crystal2DLattice.square(a=1.0)
    elements = ['H', 'He', 'B', 'C', 'N', 'O', 'Ar']
    for element in elements:
        xatom = BasisAtom(element, lattice=lattice)
        assert_equal(xatom.element, element)
        assert_equal(xatom.m, xatom.mass)

    xatom = BasisAtom()
    for c in ('x', 'y', 'z'):
        assert_equal(getattr(xatom, c), 0.0)


def test3():
    lattice = Crystal2DLattice.square(a=1.0)
    atom = BasisAtom('C', lattice=lattice)
    print(atom)


def test4():
    lattice = Crystal3DLattice.cubic(a=5.0)
    basis = BasisAtoms(atoms=['C', 'C'])
    print(basis)
    print(basis.lattice)
    assert_true(basis.lattice is None)
    basis.lattice = lattice
    assert_true(isinstance(basis.lattice, Crystal3DLattice))


def test5():
    atoms = \
        BasisAtoms(atoms=generate_atoms(elements='periodic_table').symbols)
    print(atoms[:5])

    atoms_slice = atoms[5:10]

    atoms[:5] = atoms_slice
    assert_equal(atoms[:5], atoms[5:10])
    print(atoms[:5])

    atoms[:5] = ['C', 'H', 'N', 'Ar', 'He']
    print(atoms[:8])

    atoms[0] = 'Au'
    print(atoms[:2])


if __name__ == '__main__':
    nose.runmodule()
