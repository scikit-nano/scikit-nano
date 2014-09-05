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

from ._base import GeneratorAtoms as Atoms
#from sknano.core.math import Vector

__all__ = ['NanotubeBundleGeneratorMixin']


class NanotubeBundleGeneratorMixin(object):

    def generate_bundle(self):
        self._atomsobj0 = copy.deepcopy(self.structure_atoms)
        self.structure_atoms = Atoms()
        for dr in self.bundle_coords:
            atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
            atomsobj.center_CM()
            atomsobj.translate(dr)
            self.structure_atoms.extend(atomsobj)
