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
from sknano.core.math import Point, transformation_matrix, zhat, \
    rotation_matrix
from sknano.core.refdata import aCC, element_data

r_CC_vdw = element_data['C']['VanDerWaalsRadius']


def test1():
    dlattice = Crystal2DLattice(a=4.0, b=8.0, gamma=120)
    orientation_matrix = rotation_matrix(angle=np.pi/6, axis=zhat)
    rlattice = \
        Reciprocal2DLattice(a_star=dlattice.reciprocal_lattice.a_star,
                            b_star=dlattice.reciprocal_lattice.b_star,
                            gamma_star=dlattice.reciprocal_lattice.gamma_star,
                            orientation_matrix=orientation_matrix)
    print('\ndlattice.matrix:\n{}'.format(dlattice.matrix))
    print('\nrlattice.matrix:\n{}'.format(rlattice.matrix))

    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          dlattice.reciprocal_lattice.matrix))
    print('\nrlattice.reciprocal_lattice.matrix:\n{}'.format(
          rlattice.reciprocal_lattice.matrix))

    assert_true(np.allclose(dlattice.matrix,
                            rlattice.reciprocal_lattice.matrix))
    assert_true(np.allclose(dlattice.reciprocal_lattice.matrix,
                            rlattice.matrix))


def test2():
    a = np.sqrt(3) * aCC
    latt = Crystal2DLattice(a=a, b=a, gamma=120)
    hexlatt = Crystal2DLattice.hexagonal(a)
    assert_equal(latt, hexlatt)


def test3():
    a = np.sqrt(3) * aCC
    latt = Crystal2DLattice(a=a, b=a, gamma=90)
    square = Crystal2DLattice.square(a)
    assert_equal(latt, square)


def test4():
    a = np.sqrt(3) * aCC
    latt = Crystal2DLattice(a=a, b=a, gamma=60)
    a1 = latt.a1
    a2 = latt.a2

    rotated_a1 = a1.copy()
    rotated_a2 = a2.copy()
    xfrm = transformation_matrix(angle=-np.pi / 6)
    rotated_a1.rotate(transform_matrix=xfrm)
    rotated_a2.rotate(transform_matrix=xfrm)

    latt.rotate(angle=-np.pi / 6)

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
    dlattice = Crystal2DLattice(a=a, b=a, gamma=60)
    rlattice = \
        Reciprocal2DLattice(cell_matrix=dlattice.reciprocal_lattice.matrix)
    assert_equal(dlattice, rlattice.reciprocal_lattice)
    assert_equal(dlattice.reciprocal_lattice, rlattice)


def test6():
    a = np.sqrt(3) * aCC
    l1 = Crystal2DLattice.square(a)
    l2 = Crystal2DLattice.square(2 * a)
    assert_true(l1 < l2)
    assert_true(np.allclose(2 * l1.a, l2.a))


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
    dlattice = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
                                alpha=90, beta=90, gamma=120)

    rlattice = \
        Reciprocal3DLattice(cell_matrix=dlattice.reciprocal_lattice.matrix)

    # print('\n{}'.format(np.linalg.inv(dlattice.matrix)))
    # print(np.linalg.inv(dlattice.matrix) * dlattice.matrix)
    # print(np.linalg.inv(dlattice.matrix).T)
    # print(dlattice.reciprocal_lattice.matrix)
    # print(dlattice.reciprocal_lattice.b1)
    # print(dlattice.reciprocal_lattice.b2)
    # print(dlattice.reciprocal_lattice.b3)
    # print(rlattice.matrix)
    assert_equal(dlattice, rlattice.reciprocal_lattice)
    assert_equal(dlattice.reciprocal_lattice, rlattice)
    assert_equal(dlattice.reciprocal_lattice.b1, rlattice.b1)
    assert_equal(dlattice.reciprocal_lattice.b2, rlattice.b2)
    assert_equal(dlattice.reciprocal_lattice.b3, rlattice.b3)

    assert_true(np.allclose(dlattice.reciprocal_lattice.matrix,
                            rlattice.matrix))
    assert_true(np.allclose(dlattice.reciprocal_lattice.metric_tensor,
                            rlattice.metric_tensor))

    assert_true(np.allclose(dlattice.ortho_matrix,
                            rlattice.reciprocal_lattice.ortho_matrix))
    assert_true(np.allclose(dlattice.reciprocal_lattice.ortho_matrix,
                            rlattice.ortho_matrix))
    assert_true(np.allclose(dlattice.matrix * rlattice.matrix.T, np.eye(3)))
    assert_true(np.allclose(dlattice.metric_tensor * rlattice.metric_tensor,
                            np.eye(3)))

    pmg_dlattice = pmg.Lattice(dlattice.matrix)
    print('\npmg_dlattice.matrix:\n{}'.format(
          rezero_array(pmg_dlattice.matrix)))
    print('\ndlattice.matrix:\n{}'.format(rezero_array(dlattice.matrix)))
    print('\npmg_dlattice.reciprocal_lattice.matrix:\n{}'.format(
          rezero_array(pmg_dlattice.reciprocal_lattice.matrix)))
    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          rezero_array(dlattice.reciprocal_lattice.matrix)))
    assert_true(np.allclose(dlattice.a, pmg_dlattice.a))
    assert_true(np.allclose(dlattice.b, pmg_dlattice.b))
    assert_true(np.allclose(dlattice.c, pmg_dlattice.c))
    assert_true(np.allclose(np.asarray(dlattice.a1), pmg_dlattice.matrix[0]))
    assert_true(np.allclose(np.asarray(dlattice.a2), pmg_dlattice.matrix[1]))
    assert_true(np.allclose(np.asarray(dlattice.a3), pmg_dlattice.matrix[2]))

    assert_true(np.allclose(dlattice.matrix, pmg_dlattice.matrix))

    pmg_dlattice = pmg.Lattice.from_parameters(dlattice.a,
                                               dlattice.b,
                                               dlattice.c,
                                               dlattice.alpha,
                                               dlattice.beta,
                                               dlattice.gamma)
    print('\npmg_dlattice.matrix:\n{}'.format(
          rezero_array(pmg_dlattice.matrix)))
    print('\ndlattice.matrix:\n{}'.format(rezero_array(dlattice.matrix)))
    print('\npmg_dlattice.reciprocal_lattice.matrix:\n{}'.format(
          rezero_array(pmg_dlattice.reciprocal_lattice.matrix)))
    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          rezero_array(dlattice.reciprocal_lattice.matrix)))
    assert_true(np.allclose(dlattice.a, pmg_dlattice.a))
    assert_true(np.allclose(dlattice.b, pmg_dlattice.b))
    assert_true(np.allclose(dlattice.c, pmg_dlattice.c))
    assert_true(np.allclose(np.asarray(dlattice.a1), pmg_dlattice.matrix[0]))
    assert_true(np.allclose(np.asarray(dlattice.a2), pmg_dlattice.matrix[1]))
    assert_true(np.allclose(np.asarray(dlattice.a3), pmg_dlattice.matrix[2]))

    assert_true(np.allclose(dlattice.matrix, pmg_dlattice.matrix))

    assert_true(np.allclose(dlattice.ortho_matrix.T,
                            pmg.Lattice(dlattice.ortho_matrix.T).matrix))

    assert_true(np.allclose(pmg.Lattice(dlattice.matrix).metric_tensor,
                            dlattice.metric_tensor))
    assert_true(np.allclose(pmg.Lattice(dlattice.ortho_matrix.T).metric_tensor,
                            dlattice.metric_tensor))
    assert_true(np.allclose(pmg.Lattice(dlattice.matrix).matrix,
                            dlattice.matrix))
    assert_true(np.allclose(pmg.Lattice(dlattice.ortho_matrix.T).matrix,
                            dlattice.matrix))

    assert_true(np.allclose(pmg.Lattice(dlattice.matrix).matrix,
                            rlattice.reciprocal_lattice.matrix))

    pmg_rlattice = \
        pmg.Lattice.from_parameters(dlattice.reciprocal_lattice.a_star,
                                    dlattice.reciprocal_lattice.b_star,
                                    dlattice.reciprocal_lattice.c_star,
                                    dlattice.reciprocal_lattice.alpha_star,
                                    dlattice.reciprocal_lattice.beta_star,
                                    dlattice.reciprocal_lattice.gamma_star)
    print('\npmg_rlattice:\n{}'.format(rezero_array(pmg_rlattice.matrix)))
    print('\nrlattice:\n{}'.format(rezero_array(rlattice.matrix)))
    assert_true(np.allclose(dlattice.reciprocal_lattice.matrix,
                            rlattice.matrix))
    # assert_true(np.allclose(dlattice.reciprocal_lattice.matrix,
    #                         pmg_rlattice.matrix))
    # assert_true(np.allclose(rlattice.matrix,
    #                         pmg_rlattice.matrix))


def test11():
    a = np.sqrt(3) * aCC
    dlattice = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
                                alpha=90, beta=90, gamma=120)

    orientation_matrix = rotation_matrix(angle=np.pi / 6, axis=zhat)
    rlattice = \
        Reciprocal3DLattice(a_star=dlattice.reciprocal_lattice.a_star,
                            b_star=dlattice.reciprocal_lattice.b_star,
                            c_star=dlattice.reciprocal_lattice.c_star,
                            alpha_star=dlattice.reciprocal_lattice.alpha_star,
                            beta_star=dlattice.reciprocal_lattice.beta_star,
                            gamma_star=dlattice.reciprocal_lattice.gamma_star,
                            orientation_matrix=orientation_matrix)

    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          dlattice.reciprocal_lattice.matrix))
    print('\nrlattice.matrix:\n{}'.format(rlattice.matrix))

    print('\ndlattice.matrix:\n{}'.format(dlattice.matrix))
    print('\nrlattice.reciprocal_lattice.matrix:\n{}'.format(
          rlattice.reciprocal_lattice.matrix))

    # assert_equal(dlattice.reciprocal_lattice, rlattice)
    assert_true(np.allclose(dlattice.reciprocal_lattice.matrix,
                            rlattice.matrix))
    assert_true(np.allclose(dlattice.matrix,
                            rlattice.reciprocal_lattice.matrix))

    # print('\n{}'.format(np.linalg.inv(dlattice.matrix)))
    # print(np.linalg.inv(dlattice.matrix) * dlattice.matrix)
    # print(np.linalg.inv(dlattice.matrix).T)
    # print(dlattice.reciprocal_lattice.matrix)
    # print(dlattice.reciprocal_lattice.b1)
    # print(dlattice.reciprocal_lattice.b2)
    # print(dlattice.reciprocal_lattice.b3)
    # print(rlattice.matrix)
    # print(rlattice.b1)
    # print(rlattice.b2)
    # print(rlattice.b3)

    assert_equal(dlattice.reciprocal_lattice, rlattice)
    assert_equal(dlattice.reciprocal_lattice.b1, rlattice.b1)
    assert_equal(dlattice.reciprocal_lattice.b2, rlattice.b2)
    assert_equal(dlattice.reciprocal_lattice.b3, rlattice.b3)

    assert_true(np.allclose(dlattice.reciprocal_lattice.matrix,
                            rlattice.matrix))
    assert_true(np.allclose(dlattice.reciprocal_lattice.metric_tensor,
                            rlattice.metric_tensor))

    pmg_dlattice = pmg.Lattice(dlattice.matrix)
    print('\npmg_dlattice.matrix:\n{}'.format(pmg_dlattice.matrix))
    print('\ndlattice.matrix:\n{}'.format(dlattice.matrix))
    print('\npmg_dlattice.reciprocal_lattice.matrix:\n{}'.format(
          pmg_dlattice.reciprocal_lattice.matrix))
    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          dlattice.reciprocal_lattice.matrix))
    assert_true(np.allclose(dlattice.a, pmg_dlattice.a))
    assert_true(np.allclose(dlattice.b, pmg_dlattice.b))
    assert_true(np.allclose(dlattice.c, pmg_dlattice.c))
    assert_true(np.allclose(np.asarray(dlattice.a1), pmg_dlattice.matrix[0]))
    assert_true(np.allclose(np.asarray(dlattice.a2), pmg_dlattice.matrix[1]))
    assert_true(np.allclose(np.asarray(dlattice.a3), pmg_dlattice.matrix[2]))

    assert_true(np.allclose(dlattice.matrix, pmg_dlattice.matrix))

    assert_true(np.allclose(dlattice.ortho_matrix.T,
                            pmg.Lattice(dlattice.ortho_matrix.T).matrix))

    assert_true(np.allclose(pmg.Lattice(dlattice.matrix).metric_tensor,
                            dlattice.metric_tensor))
    assert_true(np.allclose(pmg.Lattice(dlattice.ortho_matrix.T).metric_tensor,
                            dlattice.metric_tensor))
    assert_true(np.allclose(pmg.Lattice(dlattice.matrix).matrix,
                            dlattice.matrix))
    assert_true(np.allclose(pmg.Lattice(dlattice.ortho_matrix.T).matrix,
                            dlattice.matrix))

    # print('latt_from_inv_latt_matrix_transpose:\n{}\n'.format(
    #       latt_from_inv_latt_matrix_transpose.matrix))

    # print('rlattice.reciprocal_lattice.matrix:\n{}'.format(
    #       rlattice.reciprocal_lattice.matrix))

    assert_true(np.allclose(pmg.Lattice(dlattice.matrix).matrix,
                            rlattice.reciprocal_lattice.matrix))

    # print('rlattice.reciprocal_lattice.matrix:\n{}'.format(
    #       rlattice.reciprocal_lattice.matrix))
    # print('pmg.Lattice(np.linalg.inv(dlattice.matrix).T).'
    #       'reciprocal_lattice_crystallographic.matrix:\n{}'
    #       .format(pmg.Lattice(np.linalg.inv(dlattice.matrix).T)
    #               .reciprocal_lattice_crystallographic.matrix))
    # print('pmg.Lattice.from_parameters(...):\n{}\n'.format(
    #       rezero_array(pmg.Lattice.from_parameters(
    #                    dlattice.reciprocal_lattice.a_star,
    #                    dlattice.reciprocal_lattice.b_star,
    #                    dlattice.reciprocal_lattice.c_star,
    #                    dlattice.reciprocal_lattice.alpha_star,
    #                    dlattice.reciprocal_lattice.beta_star,
    #                    dlattice.reciprocal_lattice.gamma_star)
    #                    .reciprocal_lattice_crystallographic.matrix)))

    # print('dlattice.ortho_matrix:\n{}'.format(dlattice.ortho_matrix))
    # print('np.linalg.inv(dlattice.matrix).T:\n{}'.format(np.linalg.inv(
    #       dlattice.matrix).T))
    # print('rlattice.reciprocal_lattice.ortho_matrix:\n{}'
    #       .format(rlattice.reciprocal_lattice.ortho_matrix))
    # print('rlattice.reciprocal_lattice.ortho_matrix:\n{}\n'.
    #       format(rlattice.reciprocal_lattice.ortho_matrix))

    # print('dlattice.reciprocal_lattice.ortho_matrix:\n{}'.format(
    #       dlattice.reciprocal_lattice.ortho_matrix))
    # print('rlattice.ortho_matrix:\n{}'.format(
    #       rlattice.ortho_matrix))
    # print('rlattice.ortho_matrix:\n{}\n'.format(
    #       rlattice.ortho_matrix))

    assert_true(np.allclose(dlattice.ortho_matrix,
                            rlattice.reciprocal_lattice.
                            ortho_matrix))
    assert_true(np.allclose(dlattice.reciprocal_lattice.ortho_matrix,
                            rlattice.ortho_matrix))
    assert_equal(dlattice, rlattice.reciprocal_lattice)
    assert_equal(dlattice.reciprocal_lattice, rlattice)

    print(dlattice.matrix * rlattice.matrix.T)
    assert_true(np.allclose(dlattice.matrix * rlattice.matrix.T,
                np.eye(3)))


def test12():
    a = np.sqrt(3) * aCC
    cubic_latt = Crystal3DLattice(a=a, b=a, c=a, alpha=90, beta=90, gamma=90)
    assert_equal(cubic_latt, Crystal3DLattice.cubic(a))


def test13():
    latt = Crystal3DLattice(a=4.0, b=4.0, c=4.0,
                            alpha=90, beta=90, gamma=90)
    print(latt)
    p = [2.1, 0.9, 0.5]
    assert_true(np.allclose(latt.wrap_fractional_coordinate(p),
                Point((0.1, 0.9, 0.5))))


def test14():
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


def test15():
    latt = Crystal3DLattice(a=4.0, b=4.0, c=4.0,
                            alpha=90, beta=90, gamma=90)
    print(latt)

    recip_latt = Crystal3DLattice(a1=latt.b1, a2=latt.b2, a3=latt.b3)
    print(recip_latt)

    assert_equal(latt.a1, recip_latt.b1)
    assert_equal(latt.a2, recip_latt.b2)
    assert_equal(latt.a3, recip_latt.b3)


def test16():
    a = np.sqrt(3) * aCC
    assert_true(Crystal3DLattice.cubic(a) < Crystal3DLattice.cubic(2 * a))


def test17():
    dlattice = Crystal3DLattice.hexagonal(np.sqrt(3) * aCC, r_CC_vdw)
    rlattice = dlattice.reciprocal_lattice
    assert_true(np.allclose(dlattice.volume ** 2,
                            np.linalg.det(dlattice.metric_tensor)))

    print(dlattice.volume ** 2)
    print(np.linalg.det(dlattice.metric_tensor))

    print(rlattice.volume ** 2)
    print(np.linalg.det(rlattice.metric_tensor))
    assert_true(np.allclose(rlattice.volume ** 2,
                            np.linalg.det(rlattice.metric_tensor)))

    print(rlattice.metric_tensor)
    print(rlattice.metric_tensor.T)
    assert_true(np.allclose(dlattice.volume * rlattice.volume, 1.0))
    assert_true(np.allclose(rlattice.metric_tensor * dlattice.metric_tensor.T,
                            np.asmatrix(np.eye(3))))


def test18():
    a = np.sqrt(3) * aCC
    dlattice = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
                                alpha=90, beta=90, gamma=120)
    dlattice_from_matrix = Crystal3DLattice(cell_matrix=dlattice.matrix)
    assert_equal(dlattice, dlattice_from_matrix)
    assert_true(np.allclose(dlattice.matrix, dlattice_from_matrix.matrix))


def test19():
    a = np.sqrt(3) * aCC
    dlattice = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
                                alpha=90, beta=90, gamma=120)
    print('\ndlattice.matrix:\n{}'.format(dlattice.matrix))
    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          dlattice.reciprocal_lattice.matrix))

    a1 = dlattice.a1
    a2 = dlattice.a2
    a3 = dlattice.a3
    V = dlattice.volume

    b1 = a2.cross(a3) / V
    b2 = a3.cross(a1) / V
    b3 = a1.cross(a2) / V
    a_star = b1.length
    b_star = b2.length
    c_star = b3.length

    alpha_star = np.degrees(b2.angle(b3))
    beta_star = np.degrees(b3.angle(b1))
    gamma_star = np.degrees(b1.angle(b2))

    assert_true(np.allclose(a_star, dlattice.reciprocal_lattice.a_star))
    assert_true(np.allclose(b_star, dlattice.reciprocal_lattice.b_star))
    assert_true(np.allclose(c_star, dlattice.reciprocal_lattice.c_star))
    assert_true(np.allclose(alpha_star,
                            dlattice.reciprocal_lattice.alpha_star))
    assert_true(np.allclose(beta_star,
                            dlattice.reciprocal_lattice.beta_star))
    assert_true(np.allclose(gamma_star,
                            dlattice.reciprocal_lattice.gamma_star))

    rlattice = Reciprocal3DLattice(a_star=a_star, b_star=b_star, c_star=c_star,
                                   alpha_star=alpha_star, beta_star=beta_star,
                                   gamma_star=gamma_star)
    print('\nrlattice.matrix:\n{}'.format(rlattice.matrix))

    print('\nrlattice.reciprocal_lattice.matrix:\n{}'.format(
          rlattice.reciprocal_lattice.matrix))

    rlattice_from_rlattice_matrix = \
        Reciprocal3DLattice(cell_matrix=rlattice.matrix)
    print('\nrlattice_from_rlattice_matrix.matrix:\n{}'.format(
          rlattice_from_rlattice_matrix.matrix))
    print('\nrlattice_from_rlattice_matrix.reciprocal_lattice.matrix:\n{}'.
          format(rlattice_from_rlattice_matrix.reciprocal_lattice.matrix))

    assert_equal(rlattice, rlattice_from_rlattice_matrix)
    assert_true(np.allclose(rlattice.matrix,
                            rlattice_from_rlattice_matrix.matrix))

    rlattice_from_dlattice_rlattice_matrix = \
        Reciprocal3DLattice(cell_matrix=dlattice.reciprocal_lattice.matrix)
    print('\nrlattice_from_dlattice_rlattice_matrix.matrix:\n{}'.format(
          rlattice_from_dlattice_rlattice_matrix.matrix))
    print('\nrlattice_from_dlattice_rlattice_matrix.reciprocal_lattice.'
          'matrix:\n{}'.format(rlattice_from_dlattice_rlattice_matrix.
                               reciprocal_lattice.matrix))


def test20():
    a = np.sqrt(3) * aCC
    orientation_matrix = rotation_matrix(angle=-np.pi / 6, axis=zhat)
    dlattice = Crystal3DLattice(a=a, b=a, c=2 * r_CC_vdw,
                                alpha=90, beta=90, gamma=120,
                                orientation_matrix=orientation_matrix)
    rlattice = \
        Reciprocal3DLattice(a_star=dlattice.reciprocal_lattice.a_star,
                            b_star=dlattice.reciprocal_lattice.b_star,
                            c_star=dlattice.reciprocal_lattice.c_star,
                            alpha_star=dlattice.reciprocal_lattice.alpha_star,
                            beta_star=dlattice.reciprocal_lattice.beta_star,
                            gamma_star=dlattice.reciprocal_lattice.gamma_star)

    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          dlattice.reciprocal_lattice.matrix))
    print('\nrlattice.matrix:\n{}'.format(rlattice.matrix))

    print('\ndlattice.matrix:\n{}'.format(dlattice.matrix))
    print('\nrlattice.reciprocal_lattice.matrix:\n{}'.format(
          rlattice.reciprocal_lattice.matrix))
    print('\ndlattice.reciprocal_lattice.matrix:\n{}'.format(
          dlattice.reciprocal_lattice.matrix))

    pmg_rlattice = \
        pmg.Lattice.from_parameters(dlattice.reciprocal_lattice.a_star,
                                    dlattice.reciprocal_lattice.b_star,
                                    dlattice.reciprocal_lattice.c_star,
                                    dlattice.reciprocal_lattice.alpha_star,
                                    dlattice.reciprocal_lattice.beta_star,
                                    dlattice.reciprocal_lattice.gamma_star)

    print('\npmg_rlattice.matrix:\n{}'.format(pmg_rlattice.matrix))
    print('\nrlattice.matrix:\n{}'.format(rlattice.matrix))
    print('\npmg_rlattice.reciprocal_lattice.matrix:\n{}'.format(
          pmg_rlattice.reciprocal_lattice.matrix))
    print('\nrlattice.reciprocal_lattice.matrix:\n{}'.format(
          rlattice.reciprocal_lattice.matrix))

    assert_true(np.allclose(rlattice.matrix,
                            pmg_rlattice.matrix))
    assert_true(np.allclose(dlattice.reciprocal_lattice.matrix,
                            pmg_rlattice.matrix))
    assert_true(np.allclose(rlattice.matrix,
                            pmg_rlattice.matrix))


if __name__ == '__main__':
    nose.runmodule()
