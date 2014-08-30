# -*- coding: utf-8 -*-
"""
==============================================================================
Bundle Generators (:mod:`sknano.generators._nanotube_bundle_generators`)
==============================================================================

.. currentmodule:: sknano.generators._nanotube_bundle_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy
#import itertools

from sknano.core.atoms import XAtoms as Atoms
from sknano.core.math import Vector

__all__ = ['NanotubeBundleGeneratorMixin']


class NanotubeBundleGeneratorMixin(object):

    def generate_hexagonal_bundle(self):
        nrows = max(self._nx, self._ny, 3)
        if nrows % 2 != 1:
            nrows += 1

        ntubes_per_end_rows = int((nrows + 1) / 2)

        row = 0
        ntubes_per_row = nrows
        while ntubes_per_row >= ntubes_per_end_rows:
            if row == 0:
                for n in xrange(ntubes_per_row):
                    atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                    atomsobj.center_CM()
                    dr = n * self._r1
                    atomsobj.translate(dr)
                    self._structure_atoms.extend(atomsobj)
                    self._Ntubes += 1
            else:
                for nx in xrange(ntubes_per_row):
                    for ny in (-row, row):
                        atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                        atomsobj.center_CM()
                        dr = Vector()
                        dr.x = abs(ny * self._r2.x)
                        dr.y = ny * self._r2.y
                        dr = nx * self._r1 + dr
                        atomsobj.translate(dr)
                        self._structure_atoms.extend(atomsobj)
                        self._Ntubes += 1
            row += 1
            ntubes_per_row = nrows - row

    def generate_rectangular_bundle(self):
        Lx = 10 * self._Lx
        for nx in xrange(self._nx):
            for ny in xrange(self._ny):
                atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                atomsobj.center_CM()
                dr = nx * self._r1 + ny * self._r2
                while dr.x < 0:
                    dr.x += Lx
                atomsobj.translate(dr)
                self._structure_atoms.extend(atomsobj)
                self._Ntubes += 1

    def generate_square_bundle(self):
        pass

    def generate_triangular_bundle(self):
        pass

    def generate_bundle(self):
        self._Ntubes = 0
        self._atomsobj0 = copy.deepcopy(self._structure_atoms)
        self._structure_atoms = Atoms()

        if self._bundle_geometry == 'hexagon':
            self.generate_hexagonal_bundle()
        elif self._bundle_geometry == 'rectangle':
            self.generate_rectangular_bundle()
        elif self._bundle_geometry == 'square':
            self.generate_square_bundle()
        elif self._bundle_geometry == 'triangle':
            self.generate_triangular_bundle()
        else:
            for nx in xrange(self._nx):
                for ny in xrange(self._ny):
                    atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                    atomsobj.center_CM()
                    dr = nx * self._r1 + ny * self._r2
                    atomsobj.translate(dr)
                    self._structure_atoms.extend(atomsobj)
                    self._Ntubes += 1
