#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_true, assert_equal
import numpy as np

from sknano.core.geometric_regions import Parallelepiped, Cuboid, \
    generate_bounding_box
from sknano.core.math import Point, Vector
from sknano.generators import GrapheneGenerator


def test1():
    print('generating graphene structure')
    graphene = GrapheneGenerator(armchair_edge_length=5,
                                 zigzag_edge_length=5)
    lattice = graphene.lattice
    print('graphene.bounds:\n{}'.format(graphene.bounds))
    print('graphene.centroid:\n{}'.format(graphene.centroid))
    print('graphene.lattice:\n{}'.format(lattice))
    print('graphene.lattice.a1:\n{}'.format(lattice.a1))
    print('graphene.lattice.a2:\n{}'.format(lattice.a2))
    print('graphene.lattice.a3:\n{}'.format(lattice.a3))
    print('graphene.lattice.orientation_matrix:\n{}'.format(
          lattice.orientation_matrix))
    print('rotating graphene')
    graphene.rotate(angle=-np.pi/2, axis='x')
    print('graphene.bounds:\n{}'.format(graphene.bounds))
    print('graphene.centroid:\n{}'.format(graphene.centroid))
    print('graphene.lattice:\n{}'.format(lattice))
    print('graphene.lattice.a1:\n{}'.format(lattice.a1))
    print('graphene.lattice.a2:\n{}'.format(lattice.a2))
    print('graphene.lattice.a3:\n{}'.format(lattice.a3))
    print('graphene.lattice.orientation_matrix:\n{}'.format(
          lattice.orientation_matrix))

    assert_true(np.allclose(lattice.angles, 3 * [90.0]))
    lattice_region = Cuboid(pmax=lattice.lengths)

    # lattice_region = Parallelepiped(u=lattice.a * xhat,
    #                                 v=lattice.b * yhat,
    #                                 w=lattice.c * zhat)
    assert_equal(lattice_region.a, lattice.a)
    assert_equal(lattice_region.b, lattice.b)
    assert_equal(lattice_region.c, lattice.c)
    print('lattice_region:\n{}'.format(lattice_region))
    print('lattice_region.centroid:\n{}'.format(lattice_region.centroid))

    print('\nrotating lattice_region')
    lattice_region.rotate(transform_matrix=lattice.orientation_matrix)
    # assert_equal(lattice_region.a, lattice.a)
    # assert_equal(lattice_region.b, lattice.b)
    # assert_equal(lattice_region.c, lattice.c)
    print('lattice_region:\n{}'.format(lattice_region))
    print('lattice_region.centroid:\n{}'.format(lattice_region.centroid))

    print('\ncentering lattice_region on graphene centroid')
    tvec = Vector(Point(graphene.centroid) - lattice_region.centroid)
    lattice_region.translate(tvec)
    # assert_equal(lattice_region.a, lattice.a)
    # assert_equal(lattice_region.b, lattice.b)
    # assert_equal(lattice_region.c, lattice.c)
    print('lattice_region:\n{}'.format(lattice_region))
    print('lattice_region.centroid:\n{}'.format(lattice_region.centroid))

    bounding_box = generate_bounding_box(from_lattice=lattice,
                                         center=graphene.centroid,
                                         verbose=True)
    print('bounding_box:\n{}'.format(bounding_box))
    assert_equal(bounding_box, lattice_region)
    print('lattice_region.lengths: {}, {}, {}'.format(
          lattice_region.a, lattice_region.b, lattice_region.c))


def test2():
    graphene = GrapheneGenerator(armchair_edge_length=5,
                                 zigzag_edge_length=5)
    graphene.rotate(angle=-np.pi/2, axis='x')
    lattice = graphene.lattice

    bounding_box1 = Cuboid()
    bounding_box2 = Cuboid()

    lattice_region1 = Cuboid(pmax=lattice.lengths)
    bounding_box1.pmin = lattice_region1.pmin
    bounding_box1.pmax = lattice_region1.pmax

    a, b, c = lattice.lengths
    cos_alpha, cos_beta, cos_gamma = np.cos(np.radians(lattice.angles))
    lx = a
    xy = b * cos_gamma
    xz = c * cos_beta
    ly = np.sqrt(b ** 2 - xy ** 2)
    yz = (b * c * cos_alpha - xy * xz) / ly
    lz = np.sqrt(c ** 2 - xz ** 2 - yz ** 2)

    lattice_region2 = \
        Parallelepiped(u=Vector(lattice.ortho_matrix[:, 0].A.flatten()),
                       v=Vector(lattice.ortho_matrix[:, 1].A.flatten()),
                       w=Vector(lattice.ortho_matrix[:, 2].A.flatten()))

    xlo, ylo, zlo = lattice_region2.o
    print('xy={}, xz={}, yz={}'.format(xy, xz, yz))
    print('lx={}, ly={}, lz={}'.format(lx, ly, lz))
    print('xlo={}, ylo={}, zlo={}'.format(xlo, ylo, zlo))
    xlo_bound = xlo + min(0.0, xy, xz, xy + xz)
    xhi_bound = xlo + lx + max(0.0, xy, xz, xy + xz)
    ylo_bound = ylo + min(0.0, yz)
    yhi_bound = ylo + ly + max(0.0, yz)
    zlo_bound = zlo
    zhi_bound = zlo + lz
    bounding_box2.pmin = [xlo_bound, ylo_bound, zlo_bound]
    bounding_box2.pmax = [xhi_bound, yhi_bound, zhi_bound]

    print(bounding_box1)
    print(bounding_box2)

    [bounding_box.rotate(transform_matrix=lattice.orientation_matrix)
     for bounding_box in (bounding_box1, bounding_box2)]

    [bounding_box.translate(Vector(Point(graphene.centroid) -
                                   bounding_box.centroid))
     for bounding_box in (bounding_box1, bounding_box2)]

    [assert_true(bounding_box.pmin <= bounding_box.pmax) for bounding_box
     in (bounding_box1, bounding_box2)]

    assert_equal(bounding_box1, bounding_box2)

    print(bounding_box1)
    print(bounding_box2)

if __name__ == '__main__':
    nose.runmodule()
