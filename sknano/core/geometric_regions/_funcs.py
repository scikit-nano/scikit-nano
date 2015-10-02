# -*- coding: utf-8 -*-
"""
===============================================================================
Helper functions (:mod:`sknano.core.geometric_regions._funcs`)
===============================================================================

.. currentmodule:: sknano.core.geometric_regions._funcs

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import numpy as np

__all__ = ['generate_bounding_box']

from sknano.core.math import Point, Vector
from ._3D_regions import Cuboid, Parallelepiped


def generate_bounding_box(from_region=None, from_lattice=None,
                          from_array=None, center=None, verbose=False):

    if all([obj is None for obj in (from_region, from_lattice, from_array)]):
        return None

    bounding_box = Cuboid()

    if from_region is not None:
        pass
    elif from_lattice is not None:
        lattice = from_lattice
        if np.allclose(np.radians(lattice.angles), np.pi / 2 * np.ones(3)):
            lattice_region = Cuboid(pmax=lattice.lengths)
            bounding_box.pmin = lattice_region.pmin
            bounding_box.pmax = lattice_region.pmax
        else:
            a, b, c = lattice.lengths
            cos_alpha, cos_beta, cos_gamma = np.cos(np.radians(lattice.angles))
            lx = a
            xy = b * cos_gamma
            xz = c * cos_beta
            ly = np.sqrt(b ** 2 - xy ** 2)
            yz = (b * c * cos_alpha - xy * xz) / ly
            lz = np.sqrt(c ** 2 - xz ** 2 - yz ** 2)

            lattice_region = \
                Parallelepiped(
                    u=Vector(lattice.ortho_matrix[:, 0].A.flatten()),
                    v=Vector(lattice.ortho_matrix[:, 1].A.flatten()),
                    w=Vector(lattice.ortho_matrix[:, 2].A.flatten()))

            xlo, ylo, zlo = lattice_region.o
            if verbose:
                print('xy={}, xz={}, yz={}'.format(xy, xz, yz))
                print('lx={}, ly={}, lz={}'.format(lx, ly, lz))
                print('xlo={}, ylo={}, zlo={}'.format(xlo, ylo, zlo))
            xlo_bound = xlo + min(0.0, xy, xz, xy + xz)
            xhi_bound = xlo + lx + max(0.0, xy, xz, xy + xz)
            ylo_bound = ylo + min(0.0, yz)
            yhi_bound = ylo + ly + max(0.0, yz)
            zlo_bound = zlo
            zhi_bound = zlo + lz
            bounding_box.pmin = [xlo_bound, ylo_bound, zlo_bound]
            bounding_box.pmax = [xhi_bound, yhi_bound, zhi_bound]

        if verbose:
            print(bounding_box)

        bounding_box.rotate(transform_matrix=lattice.orientation_matrix)

        if center is not None:
            tvec = Vector(Point(center) - bounding_box.centroid)
            bounding_box.translate(tvec)
        assert bounding_box.pmin <= bounding_box.pmax
    else:
        array = from_array
        for i, dim in enumerate(('x', 'y', 'z')):
            setattr(bounding_box, dim + 'min', array[:, i].min())
            setattr(bounding_box, dim + 'max', array[:, i].max())

    return bounding_box
