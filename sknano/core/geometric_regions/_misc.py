# -*- coding: utf-8 -*-
"""
===============================================================================
Helper functions (:mod:`sknano.core.geometric_regions._misc`)
===============================================================================

.. currentmodule:: sknano.core.geometric_regions._misc

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import numpy as np

__all__ = ['generate_bounding_box']

from ._3D_regions import Cuboid


def generate_bounding_box(from_region=None, from_lattice=None, from_array=None,
                          verbose=False):
    """Return a :class:`~sknano.core.geometric_regions.Cuboid` \
        representing an axis-aligned bounding box.

    Parameters
    ----------
    from_region : :class:`~sknano.core.geometric_regions.Geometric3DRegion`
    from_lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`
    from_array : :class:`~numpy:numpy.ndarray`
    verbose : :class:`~python:bool`

    Returns
    -------
    bounding_box : :class:`~sknano.core.geometric_regions.Cuboid`

    """
    if all([obj is None for obj in (from_region, from_lattice, from_array)]):
        return None

    bounding_box = Cuboid()

    if from_region is not None:
        region = from_region
        bounding_box.pmin = region.pmin
        bounding_box.pmax = region.pmax
    elif from_lattice is not None:
        lattice = from_lattice
        # if np.allclose(np.radians(lattice.angles), np.pi / 2 * np.ones(3)):
        #     lattice_region =
        #     bounding_box.pmin = lattice_region.pmin
        #     bounding_box.pmax = lattice_region.pmax
        # else:
        a, b, c = lattice.lengths
        cos_alpha, cos_beta, cos_gamma = np.cos(np.radians(lattice.angles))
        lx = a
        xy = b * cos_gamma
        xz = c * cos_beta
        ly = np.sqrt(b ** 2 - xy ** 2)
        yz = (b * c * cos_alpha - xy * xz) / ly
        lz = np.sqrt(c ** 2 - xz ** 2 - yz ** 2)

        lattice_region = lattice.region

        xlo, ylo, zlo = lattice_region.o
        if verbose:
            print('xy={}, xz={}, yz={}'.format(xy, xz, yz))
            print('lx={}, ly={}, lz={}'.format(lx, ly, lz))
            print('xlo={}, ylo={}, zlo={}'.format(xlo, ylo, zlo))
            print('lattice.region.centroid: {}'.format(
                lattice_region.centroid))
        xlo_bound = xlo + min(0.0, xy, xz, xy + xz)
        xhi_bound = xlo + lx + max(0.0, xy, xz, xy + xz)
        ylo_bound = ylo + min(0.0, yz)
        yhi_bound = ylo + ly + max(0.0, yz)
        zlo_bound = zlo
        zhi_bound = zlo + lz
        bounding_box.pmin = [xlo_bound, ylo_bound, zlo_bound]
        bounding_box.pmax = [xhi_bound, yhi_bound, zhi_bound]

        if verbose:
            print('bounding_box: {}'.format(bounding_box))
            print('orientation_matrix: {}'.format(lattice.orientation_matrix))

        assert bounding_box.pmin <= bounding_box.pmax
    else:
        array = np.asarray(from_array)
        for i, dim in enumerate(('x', 'y', 'z')):
            setattr(bounding_box, dim + 'min', array[:, i].min())
            setattr(bounding_box, dim + 'max', array[:, i].max())

    return bounding_box
