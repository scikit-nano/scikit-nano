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
                          from_array=None, center=None):

    if all([obj is None for obj in (from_region, from_lattice, from_array)]):
        return None

    bounding_box = Cuboid()

    if from_region is not None:
        pass
    elif from_lattice is not None:
        lattice = from_lattice
        a, b, c = lattice.a, lattice.b, lattice.c
        cos_alpha, cos_beta, cos_gamma = \
            lattice.cos_alpha, lattice.cos_beta, lattice.cos_gamma
        lx = a
        xy = b * cos_gamma
        xz = c * cos_beta
        ly = np.sqrt(b ** 2 - xy ** 2)
        yz = (b * c * cos_alpha - xy * xz) / ly
        lz = np.sqrt(c ** 2 - xz ** 2 - yz ** 2)

        lattice_region = \
            Parallelepiped(u=lattice.a1, v=lattice.a2, w=lattice.a3)
        if center is not None:
            tvec = Vector(Point(center) - lattice_region.centroid)
            lattice_region.translate(tvec)

        xlo, ylo, zlo = lattice_region.o
        xlo_bound = xlo + min(0.0, xy, xz, xy + xz)
        xhi_bound = xlo + lx + max(0.0, xy, xz, xy + xz)
        ylo_bound = ylo + min(0.0, yz)
        yhi_bound = ylo + ly + max(0.0, yz)
        zlo_bound = zlo
        zhi_bound = zlo + lz
        bounding_box.pmin = [xlo_bound, ylo_bound, zlo_bound]
        bounding_box.pmax = [xhi_bound, yhi_bound, zhi_bound]
    else:
        array = from_array
        for i, dim in enumerate(('x', 'y', 'z')):
            setattr(bounding_box, dim + 'min', array[:, i].min())
            setattr(bounding_box, dim + 'max', array[:, i].max())

    return bounding_box
