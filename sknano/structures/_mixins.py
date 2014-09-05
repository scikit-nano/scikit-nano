# -*- coding: utf-8 -*-
"""
==============================================================================
Mixin structure classes (:mod:`sknano.structures._mixins`)
==============================================================================

.. currentmodule:: sknano.structures._mixins

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers

import numpy as np

from sknano.core.math import Vector

__all__ = ['MWNTMixin', 'NanotubeBundleMixin', 'UnrolledSWNTMixin']


class MWNTMixin(object):
    pass


class NanotubeBundleMixin(object):

    def compute_bundle_params(self):
        """Compute/update nanotube bundle parameters."""

        self.r1.x = self.dt + self.vdw_spacing
        if self.bundle_packing is None and \
                self.bundle_geometry in ('square', 'rectangle'):
            self._bundle_packing = 'ccp'
        elif self.bundle_packing is None and \
                self.bundle_geometry in ('triangle', 'hexagon'):
            self._bundle_packing = 'hcp'

        if self.bundle_packing in ('cubic', 'ccp'):
            self.r2.y = self.r1.x
        else:
            self.r2.x = self.r1.x * np.cos(2 * np.pi / 3)
            self.r2.y = self.r1.x * np.sin(2 * np.pi / 3)
            if self.bundle_packing is None:
                self._bundle_packing = 'hcp'

        if self.bundle_geometry == 'hexagon':
            nrows = max(self.nx, self.ny, 3)
            if nrows % 2 != 1:
                nrows += 1

            ntubes_per_end_rows = int((nrows + 1) / 2)

            row = 0
            ntubes_per_row = nrows
            while ntubes_per_row >= ntubes_per_end_rows:
                if row == 0:
                    for n in xrange(ntubes_per_row):
                        dr = n * self.r1
                        self.bundle_coords.append(dr)
                else:
                    for nx in xrange(ntubes_per_row):
                        for ny in (-row, row):
                            dr = Vector()
                            dr.x = abs(ny * self.r2.x)
                            dr.y = ny * self.r2.y
                            dr = nx * self.r1 + dr
                            self.bundle_coords.append(dr)
                row += 1
                ntubes_per_row = nrows - row

        elif self.bundle_geometry == 'rectangle':
            Lx = 10 * self.Lx
            for nx in xrange(self.nx):
                for ny in xrange(self.ny):
                    dr = nx * self.r1 + ny * self.r2
                    while dr.x < 0:
                        dr.x += Lx
                    self.bundle_coords.append(dr)

        elif self.bundle_geometry == 'square':
            pass
        elif self.bundle_geometry == 'triangle':
            pass
        else:
            for nx in xrange(self.nx):
                for ny in xrange(self.ny):
                    dr = nx * self.r1 + ny * self.r2
                    self.bundle_coords.append(dr)

    @property
    def nx(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return self._nx

    @nx.setter
    def nx(self, value):
        """Set :math:`n_x`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a positive integer.')
        self._nx = int(value)

    @nx.deleter
    def nx(self):
        del self._nx

    @property
    def ny(self):
        """Number of nanotubes along the :math:`y`-axis."""
        return self._ny

    @ny.setter
    def ny(self, value):
        """Set :math:`n_y`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a positive integer.')
        self._ny = int(value)

    @ny.deleter
    def ny(self):
        del self._ny

    @property
    def Lx(self):
        return self.nx * (self.dt + self.vdw_spacing) / 10

    @property
    def Ly(self):
        return self.ny * (self.dt + self.vdw_spacing) / 10

    @property
    def bundle_packing(self):
        return self._bundle_packing

    @bundle_packing.setter
    def bundle_packing(self, value):
        if value not in ('ccp', 'hcp'):
            raise ValueError('Expected value to be `hcp` or `ccp`')
        self._bundle_packing = value
        self.compute_bundle_params()

    @bundle_packing.deleter
    def bundle_packing(self):
        del self._bundle_packing


class UnrolledSWNTMixin(object):
    @property
    def nx(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return self._nx

    @nx.setter
    def nx(self, value):
        """Set :math:`n_x`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a real positive number.')
        self._nx = int(value)

    @property
    def ny(self):
        """Number of nanotubes along the :math:`y`-axis."""
        return self._ny

    @ny.setter
    def ny(self, value):
        """Set :math:`n_y`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a real positive number.')
        self._ny = int(value)
