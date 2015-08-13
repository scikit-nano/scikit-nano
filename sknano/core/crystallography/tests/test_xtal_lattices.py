#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_equal, assert_true
import numpy as np

import pymatgen as pmg

from sknano.core import rezero_array
from sknano.core.crystallography import Crystal2DLattice, Crystal3DLattice, \
    Reciprocal2DLattice, Reciprocal3DLattice
# from sknano.core.atoms import Atom, Atoms, XAtom, XAtoms
from sknano.core.math import Point, transformation_matrix
from sknano.core.refdata import aCC, element_data

r_CC_vdw = element_data['C']['VanDerWaalsRadius']


def test1():
    latt = Crystal2DLattice(a=4.0, b=8.0, gamma=120)
    print(latt)


def test2():
    a = np.sqrt(3) * 1.42
    latt = Crystal2DLattice(a=a, b=a, gamma=120)
    hexlatt = Crystal2DLattice.hexagonal(a)
    assert_equal(latt, hexlatt)


def test3():
    a = np.sqrt(3) * 1.42
    latt = Crystal2DLattice(a=a, b=a, gamma=90)
    square = Crystal2DLattice.square(a)
    assert_equal(latt, square)


def test4():
    a = np.sqrt(3) * 1.42
    latt = Crystal2DLattice(a=a, b=a, gamma=60)
    a1 = latt.a1
    a2 = latt.a2

    rotated_a1 = a1.copy()
    rotated_a2 = a2.copy()
    xfrm = transformation_matrix(angle=-np.pi/6)
    rotated_a1.rotate(transform_matrix=xfrm)
    rotated_a2.rotate(transform_matrix=xfrm)

    latt.rotate(angle=-np.pi/6)

    assert_equal(latt.a1, rotated_a1)
    assert_equal(latt.a2, rotated_a2)
    assert_true(np.allclose(latt.orientation_matrix, xfrm))

    rotated_latt = Crystal2DLattice(a1=rotated_a1, a2=rotated_a2)
    assert_equal(rotated_a1, rotated_latt.a1)
    assert_equal(rotated_a2, rotated_latt.a2)
    assert_true(np.allclose(latt.orientation_matrix,
                            rotated_latt.orientation_matrix))


def test5():
    a = np.sqrt(3) * aCC
    latt = Crystal2DLattice(a=a, b=a, gamma=60)
    recip_latt = \
        Reciprocal2DLattice(a_star=latt.reciprocal_lattice.a_star,
                            b_star=latt.reciprocal_lattice.b_star,
                            gamma_star=latt.reciprocal_lattice.gamma_star)
    assert_equal(latt, recip_latt.reciprocal_lattice)
    assert_equal(latt.reciprocal_lattice, recip_latt)


def test6():
    a = np.sqrt(3) * aCC
    l1 = Crystal2DLattice.square(a)
    l2 = Crystal2DLattice.square(2*a)
    assert_true(l1 < l2)


def test7():
    latt = Crystal3DLattice(a=4.0, b=8.0, c=2.0, alpha=90,
                            beta=90, gamma=120)
    assert_true(np.allclose(latt.a, 4.0))
    assert_true(np.allclose(latt.b, 8.0))
    assert_true(np.allclose(latt.c, 2.0))
    assert_true(np.allclose(latt.alpha, 90.))
    assert_true(np.allclose(latt.beta, 90.))
    assert_true(np.allclose(latt.gamma, 120.))


def test8():
    a = np.sqrt(3) * aCC
    latt = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
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


def test9():
    a = np.sqrt(3) * aCC
    latt = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
                            alpha=90, beta=90, gamma=120)
    print(latt)
    hex_latt = \
        Crystal3DLattice.hexagonal(a, 2 * r_CC_vdw)
    print(hex_latt)
    assert_equal(latt, hex_latt)


def test10():
    a = np.sqrt(3) * aCC
    latt = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
                            alpha=90, beta=90, gamma=120)

    # latt_from_inv_latt_matrix_transpose = \
    #     Crystal3DLattice(cell_matrix=np.linalg.inv(latt.matrix).T)

    recip_latt_from_matrix = \
        Reciprocal3DLattice(cell_matrix=np.linalg.inv(latt.ortho_matrix).T)
    assert_equal(latt.reciprocal_lattice, recip_latt_from_matrix)
    assert_true(np.allclose(latt.reciprocal_lattice.matrix,
                            recip_latt_from_matrix.matrix))

    recip_latt_from_params = \
        Reciprocal3DLattice(a_star=latt.reciprocal_lattice.a_star,
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


def test11():
    a = np.sqrt(3) * aCC
    cubic_latt = Crystal3DLattice(a=a, b=a, c=a, alpha=90, beta=90, gamma=90)
    assert_equal(cubic_latt, Crystal3DLattice.cubic(a))


def test12():
    latt = Crystal3DLattice(a=4.0, b=4.0, c=4.0,
                            alpha=90, beta=90, gamma=90)
    print(latt)
    p = [2.1, 0.9, 0.5]
    assert_true(np.allclose(latt.wrap_fractional_coordinate(p),
                Point((0.1, 0.9, 0.5))))


def test13():
    latt = Crystal3DLattice(a=4.0, b=4.0, c=4.0,
                            alpha=90, beta=90, gamma=90)
    print(latt)

    a = latt.a1
    b = latt.a2
    c = latt.a3
    G = np.matrix([[a.dot(a), a.dot(b), a.dot(c)],
                   [b.dot(a), b.dot(b), b.dot(c)],
                   [c.dot(a), c.dot(b), c.dot(c)]])

    assert_true(np.allclose(latt.metric_tensor, G))


def test14():
    latt = Crystal3DLattice(a=4.0, b=4.0, c=4.0,
                            alpha=90, beta=90, gamma=90)
    print(latt)

    recip_latt = Crystal3DLattice(a1=latt.b1, a2=latt.b2, a3=latt.b3)
    print(recip_latt)

    assert_equal(latt.a1, recip_latt.b1)
    assert_equal(latt.a2, recip_latt.b2)
    assert_equal(latt.a3, recip_latt.b3)


def test15():
    a = np.sqrt(3) * aCC
    assert_true(Crystal3DLattice.cubic(a) < Crystal3DLattice.cubic(2 * a))


if __name__ == '__main__':
    nose.runmodule()
