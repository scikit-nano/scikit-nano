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

__all__ = ['Domain', 'generate_bounding_box']

from sknano.core import BaseClass, TabulateMixin
# from sknano.core.math import Point, Vector
from ._3D_regions import Cuboid


class Domain(TabulateMixin, BaseClass):
    """Container class for molecular dynamics simulation box metadata.

    Attributes
    ----------
    bounding_box : :class:`Cuboid`
    triclinic : :class:`~python:bool`
    xy, xz, yz : :class:`~python:float`

    """
    def __init__(self):
        self.bounding_box = Cuboid()
        self.lattice = None
        self.triclinic = False
        self.xy = self.xz = self.yz = 0.0

    def __str__(self):
        strrep = self._table_title_str()
        lattice = self.lattice
        if lattice is not None:
            strrep = '\n'.join((strrep, str(lattice)))

        # lattice_params = \
        #     self._tabulate(list(zip(('a', 'b', 'c', 'α', 'β', 'γ'),
        #                             self.lengths_and_angles)))

        return strrep

    @property
    def xlo(self):
        """Alias for :attr:`Domain.bounding_box.xmin`"""
        return self.bounding_box.xmin

    @xlo.setter
    def xlo(self, value):
        self.bounding_box.xmin = value

    @property
    def xhi(self):
        """Alias for :attr:`Domain.bounding_box.xmax`"""
        return self.bounding_box.xmax

    @xhi.setter
    def xhi(self, value):
        self.bounding_box.xmax = value

    @property
    def ylo(self):
        """Alias for :attr:`Domain.bounding_box.ymin`"""
        return self.bounding_box.ymin

    @ylo.setter
    def ylo(self, value):
        self.bounding_box.ymin = value

    @property
    def yhi(self):
        """Alias for :attr:`Domain.bounding_box.ymax`"""
        return self.bounding_box.ymax

    @yhi.setter
    def yhi(self, value):
        self.bounding_box.ymax = value

    @property
    def zlo(self):
        """Alias for :attr:`Domain.bounding_box.zmin`"""
        return self.bounding_box.zmin

    @zlo.setter
    def zlo(self, value):
        self.bounding_box.zmin = value

    @property
    def zhi(self):
        """Alias for :attr:`Domain.bounding_box.zmax`"""
        return self.bounding_box.zmax

    @zhi.setter
    def zhi(self, value):
        self.bounding_box.zmax = value

    @property
    def lx(self):
        """Alias for :attr:`Domain.bounding_box.a`."""
        return self.bounding_box.a

    @property
    def ly(self):
        """Alias for :attr:`Domain.bounding_box.a`."""
        return self.bounding_box.b

    @property
    def lz(self):
        """Alias for :attr:`Domain.bounding_box.a`."""
        return self.bounding_box.c

    @property
    def xlo_bound(self):
        """Triclinic bounding box minimum extent in the x-dimension"""
        xlo = self.xlo
        xy = self.xy
        xz = self.xz
        return xlo + min((0.0, xy, xz, xy + xz))

    @xlo_bound.setter
    def xlo_bound(self, value):
        xy = self.xy
        xz = self.xz
        self.xlo = value - min((0.0, xy, xz, xy + xz))

    @property
    def xhi_bound(self):
        """Triclinic bounding box maximum extent in the x-dimension"""
        xhi = self.xhi
        xy = self.xy
        xz = self.xz
        return xhi + max((0.0, xy, xz, xy + xz))

    @xhi_bound.setter
    def xhi_bound(self, value):
        xy = self.xy
        xz = self.xz
        self.xhi = value - max((0.0, xy, xz, xy + xz))

    @property
    def ylo_bound(self):
        """Triclinic bounding box minimum extent in the y-dimension"""
        return self.ylo + min((0.0, self.yz))

    @ylo_bound.setter
    def ylo_bound(self, value):
        self.ylo = value - min((0.0, self.yz))

    @property
    def yhi_bound(self):
        """Triclinic bounding box maximum extent in the y-dimension"""
        return self.yhi + max((0.0, self.yz))

    @yhi_bound.setter
    def yhi_bound(self, value):
        self.yhi = value - max((0.0, self.yz))

    @property
    def zlo_bound(self):
        """Triclinic bounding box minimum extent in the z-dimension"""
        return self.zlo

    @zlo_bound.setter
    def zlo_bound(self, value):
        self.zlo = value

    @property
    def zhi_bound(self):
        """Triclinic bounding box maximum extent in the z-dimension"""
        return self.zhi

    @zhi_bound.setter
    def zhi_bound(self, value):
        self.zhi = value

    def update(self, from_lattice=None, from_region=None, from_array=None,
               allow_triclinic_box=False, pad_box=False,
               pad_tol=0.01, xpad=10., ypad=10., zpad=10., verbose=False,
               **kwargs):
        """Update simulation domain attributes from lattice."""
        bounding_box = \
            generate_bounding_box(from_lattice=from_lattice,
                                  from_region=from_region,
                                  from_array=from_array,
                                  verbose=verbose)

        if pad_box and from_array is not None:
            coords = from_array
            boxpad = {'x': xpad, 'y': ypad, 'z': zpad}
            # for dim, pad in boxpad.items():
            for i, dim in enumerate(('x', 'y', 'z')):
                pad = boxpad[dim]
                dmin = dim + 'min'
                dmax = dim + 'max'
                if abs(getattr(bounding_box, dmin) -
                       coords[:, i].min()) < pad - pad_tol:
                    setattr(bounding_box, dmin,
                            getattr(bounding_box, dmin) - pad)
                if abs(getattr(bounding_box, dmax) -
                       coords[:, i].max()) < pad - pad_tol:
                    setattr(bounding_box, dmax,
                            getattr(bounding_box, dmax) + pad)

        if allow_triclinic_box and from_lattice is not None:
            self.lattice = lattice = from_lattice
            if not np.allclose(np.radians(lattice.angles),
                               np.pi / 2 * np.ones(3)):
                self.triclinic = True
                a, b, c = lattice.lengths
                cos_alpha, cos_beta, cos_gamma = \
                    np.cos(np.radians(lattice.angles))
                self.xy = xy = b * cos_gamma
                self.xz = xz = c * cos_beta
                self.yz = \
                    (b * c * cos_alpha - xy * xz) / np.sqrt(b ** 2 - xy ** 2)

                self.xlo_bound = bounding_box.xmin
                self.xhi_bound = bounding_box.xmax
                self.ylo_bound = bounding_box.ymin
                self.yhi_bound = bounding_box.ymax
                self.zlo_bound = bounding_box.zmin
                self.zhi_bound = bounding_box.zmax
        else:
            self.bounding_box = bounding_box

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict()


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
