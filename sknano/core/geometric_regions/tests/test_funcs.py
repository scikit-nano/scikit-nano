#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import assert_true, assert_equal
import numpy as np

from sknano.core.geometric_regions import Parallelepiped, Cuboid, \
    generate_bounding_box
from sknano.core.math import Point, Vector, xhat, yhat, zhat
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


if __name__ == '__main__':
    nose.runmodule()
