#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_true
import numpy as np

import pymatgen as pmg

from sknano.core import rezero_array
from sknano.core.crystallography import CrystalLattice, ReciprocalLattice
# from sknano.core.atoms import Atom, Atoms, XAtom, XAtoms
from sknano.core.math import Point, transformation_matrix
from sknano.core.refdata import aCC, element_data

r_CC_vdw = element_data['C']['VanDerWaalsRadius']


def test1():
    latt = CrystalLattice(a=4.0, b=8.0, c=2.0, alpha=90,
                          beta=90, gamma=120)
    assert_true(np.allclose(latt.a, 4.0))
    assert_true(np.allclose(latt.b, 8.0))
    assert_true(np.allclose(latt.c, 2.0))
    assert_true(np.allclose(latt.alpha, 90.))
    assert_true(np.allclose(latt.beta, 90.))
    assert_true(np.allclose(latt.gamma, 120.))


def test2():
    a = np.sqrt(3) * aCC
    latt = CrystalLattice(a=a, b=a, c=2 * r_CC_vdw,
                          alpha=90, beta=90, gamma=120)
    print(latt)
    a1 = latt.a1
    a2 = latt.a2
    a3 = latt.a3

    xfrm = transformation_matrix(angle=-np.pi / 6)

    rotangle = -np.pi / 6
    for v in (a1, a2, a3):
        v.rotate(angle=rotangle)

    latt.rotate(angle=rotangle, axis='z')
    print(latt)

    assert_equal(latt.a1, a1)
    assert_equal(latt.a2, a2)
    assert_equal(latt.a3, a3)

    assert_true(np.allclose(latt.orientation_matrix, xfrm))


def test3():
    a = np.sqrt(3) * aCC
    latt = CrystalLattice(a=a, b=a, c=2 * r_CC_vdw,
                          alpha=90, beta=90, gamma=120)
    print(latt)
    hex_latt = \
        CrystalLattice.hexagonal(a, 2 * r_CC_vdw)
    print(hex_latt)
    assert_equal(latt, hex_latt)


def test4():
    a = np.sqrt(3) * aCC
    latt = CrystalLattice(a=a, b=a, c=2 * r_CC_vdw,
                          alpha=90, beta=90, gamma=120)

    # latt_from_inv_latt_matrix_transpose = \
    #     CrystalLattice(cell_matrix=np.linalg.inv(latt.matrix).T)

    recip_latt_from_matrix = \
        ReciprocalLattice(cell_matrix=np.linalg.inv(latt.ortho_matrix).T)
    assert_equal(latt.reciprocal_lattice, recip_latt_from_matrix)
    assert_true(np.allclose(latt.reciprocal_lattice.matrix,
                            recip_latt_from_matrix.matrix))

    recip_latt_from_params = \
        ReciprocalLattice(a_star=latt.reciprocal_lattice.a_star,
                          b_star=latt.reciprocal_lattice.b_star,
                          c_star=latt.reciprocal_lattice.c_star,
                          alpha_star=latt.reciprocal_lattice.alpha_star,
                          beta_star=latt.reciprocal_lattice.beta_star,
                          gamma_star=latt.reciprocal_lattice.gamma_star)
    # print('latt:\n{}'.format(latt))
    # print('recip_latt_from_matrix:\n{}'.format(
    #       recip_latt_from_matrix))
    # print('recip_latt_from_params:\n{}'.format(
    #       recip_latt_from_params))

    pmg_latt = pmg.Lattice(latt.matrix)
    assert_true(np.allclose(latt.matrix, pmg_latt.matrix))

    # print('pmg.Lattice(latt.matrix).matrix:\n{}'.format(
    #       pmg.Lattice(latt.matrix).matrix))
    # print('latt.matrix:\n{}'.format(latt.matrix))

    # print('pmg.Lattice(latt.ortho_matrix.T).matrix:\n{}'.format(
    #       pmg.Lattice(latt.ortho_matrix.T).matrix))
    # print('latt.ortho_matrix.T:\n{}'.format(latt.ortho_matrix.T))

    assert_true(np.allclose(latt.ortho_matrix.T,
                            pmg.Lattice(latt.ortho_matrix.T).matrix))

    assert_true(np.allclose(pmg.Lattice(latt.matrix).metric_tensor,
                            latt.metric_tensor))
    assert_true(np.allclose(pmg.Lattice(latt.ortho_matrix.T).metric_tensor,
                            latt.metric_tensor))
    assert_true(np.allclose(pmg.Lattice(latt.matrix).matrix, latt.matrix))
    assert_true(np.allclose(pmg.Lattice(latt.ortho_matrix.T).matrix,
                latt.matrix))

    # print('latt_from_inv_latt_matrix_transpose:\n{}\n'.format(
    #       latt_from_inv_latt_matrix_transpose.matrix))

    # print('recip_latt_from_matrix.reciprocal_lattice.matrix:\n{}'.format(
    #       recip_latt_from_matrix.reciprocal_lattice.matrix))

    assert_true(np.allclose(pmg.Lattice(latt.matrix).matrix,
                            recip_latt_from_matrix.reciprocal_lattice.matrix))

    # print('recip_latt_from_params.reciprocal_lattice.matrix:\n{}'.format(
    #       recip_latt_from_params.reciprocal_lattice.matrix))
    # print('pmg.Lattice(np.linalg.inv(latt.matrix).T).'
    #       'reciprocal_lattice_crystallographic.matrix:\n{}'
    #       .format(pmg.Lattice(np.linalg.inv(latt.matrix).T)
    #               .reciprocal_lattice_crystallographic.matrix))
    # print('pmg.Lattice.from_parameters(...):\n{}\n'.format(
    #       rezero_array(pmg.Lattice.from_parameters(
    #                    latt.reciprocal_lattice.a_star,
    #                    latt.reciprocal_lattice.b_star,
    #                    latt.reciprocal_lattice.c_star,
    #                    latt.reciprocal_lattice.alpha_star,
    #                    latt.reciprocal_lattice.beta_star,
    #                    latt.reciprocal_lattice.gamma_star)
    #                    .reciprocal_lattice_crystallographic.matrix)))

    pmg_recip_latt_from_params = \
        pmg.Lattice.from_parameters(latt.reciprocal_lattice.a_star,
                                    latt.reciprocal_lattice.b_star,
                                    latt.reciprocal_lattice.c_star,
                                    latt.reciprocal_lattice.alpha_star,
                                    latt.reciprocal_lattice.beta_star,
                                    latt.reciprocal_lattice.gamma_star)
    assert_true(np.allclose(recip_latt_from_params.matrix,
                            pmg_recip_latt_from_params.matrix))

    assert_true(np.allclose(latt.reciprocal_lattice.matrix,
                            recip_latt_from_matrix.matrix))
    assert_true(np.allclose(latt.reciprocal_lattice.matrix,
                            recip_latt_from_params.matrix))
    assert_true(np.allclose(latt.reciprocal_lattice.matrix,
                            pmg_recip_latt_from_params.matrix))

    # print('latt.reciprocal_lattice.matrix:\n{}'.format(
    #       latt.reciprocal_lattice.matrix))
    # print('recip_latt_from_matrix.matrix:\n{}'.format(
    #       recip_latt_from_matrix.matrix))
    # print('recip_latt_from_params.matrix:\n{}'.format(
    #       recip_latt_from_params.matrix))
    # print('pmg_latt.reciprocal_lattice_crystallographic.matrix:\n{}\n'.format(
    #       pmg_latt.reciprocal_lattice_crystallographic.matrix))

    # print('latt.ortho_matrix:\n{}'.format(latt.ortho_matrix))
    # print('np.linalg.inv(latt.matrix).T:\n{}'.format(np.linalg.inv(
    #       latt.matrix).T))
    # print('recip_latt_from_matrix.reciprocal_lattice.ortho_matrix:\n{}'
    #       .format(recip_latt_from_matrix.reciprocal_lattice.ortho_matrix))
    # print('recip_latt_from_params.reciprocal_lattice.ortho_matrix:\n{}\n'.
    #       format(recip_latt_from_params.reciprocal_lattice.ortho_matrix))

    # print('latt.reciprocal_lattice:\n{}'.format(latt.reciprocal_lattice))
    # print('recip_latt_from_matrix.reciprocal_lattice:\n{}'.format(
    #       recip_latt_from_matrix.reciprocal_lattice))
    # print('recip_latt_from_params.reciprocal_lattice:\n{}'.format(
    #       recip_latt_from_params.reciprocal_lattice))

    # print('latt.reciprocal_lattice.ortho_matrix:\n{}'.format(
    #       latt.reciprocal_lattice.ortho_matrix))
    # print('recip_latt_from_matrix.ortho_matrix:\n{}'.format(
    #       recip_latt_from_matrix.ortho_matrix))
    # print('recip_latt_from_params.ortho_matrix:\n{}\n'.format(
    #       recip_latt_from_params.ortho_matrix))

    assert_true(np.allclose(latt.ortho_matrix,
                            recip_latt_from_matrix.reciprocal_lattice.
                            ortho_matrix))
    assert_true(np.allclose(latt.ortho_matrix,
                            recip_latt_from_params.reciprocal_lattice.
                            ortho_matrix))

    assert_true(np.allclose(latt.reciprocal_lattice.ortho_matrix,
                            recip_latt_from_matrix.ortho_matrix))
    assert_true(np.allclose(latt.reciprocal_lattice.ortho_matrix,
                            recip_latt_from_params.ortho_matrix))

    assert_equal(latt, recip_latt_from_params.reciprocal_lattice)
    assert_equal(latt.reciprocal_lattice, recip_latt_from_params)


def test5():
    a = np.sqrt(3) * aCC
    cubic_latt = CrystalLattice(a=a, b=a, c=a, alpha=90, beta=90, gamma=90)
    assert_equal(cubic_latt, CrystalLattice.cubic(a))


def test6():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)
    print(latt)
    p = [2.1, 0.9, 0.5]
    assert_true(np.allclose(latt.wrap_fractional_coordinate(p),
                Point((0.1, 0.9, 0.5))))


def test7():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)
    print(latt)

    a = latt.a1
    b = latt.a2
    c = latt.a3
    G = np.matrix([[a.dot(a), a.dot(b), a.dot(c)],
                   [b.dot(a), b.dot(b), b.dot(c)],
                   [c.dot(a), c.dot(b), c.dot(c)]])

    assert_true(np.allclose(latt.metric_tensor, G))


def test8():
    latt = CrystalLattice(a=4.0, b=4.0, c=4.0,
                          alpha=90, beta=90, gamma=90)
    print(latt)

    recip_latt = CrystalLattice(a1=latt.b1, a2=latt.b2, a3=latt.b3)
    print(recip_latt)

    assert_equal(latt.a1, recip_latt.b1)
    assert_equal(latt.a2, recip_latt.b2)
    assert_equal(latt.a3, recip_latt.b3)


def test9():
    a = np.sqrt(3) * aCC
    assert_true(CrystalLattice.cubic(a) < CrystalLattice.cubic(2 * a))


if __name__ == '__main__':
    nose.runmodule()
