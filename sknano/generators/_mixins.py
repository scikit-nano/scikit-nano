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

from ._base import Atoms
#from sknano.core.math import Vector

__all__ = ['CappedNanotubeGeneratorMixin',
           'NanotubeBundleGeneratorMixin']


class CappedNanotubeGeneratorMixin(object):

    def generate_endcaps(self):
        pass


class NanotubeBundleGeneratorMixin(object):

    def generate_bundle(self):
        self._atomsobj0 = copy.deepcopy(self.atoms)
        self.structure_data.clear()
        #self.atoms = Atoms()
        for dr in self.bundle_coords:
            atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
            atomsobj.center_CM()
            atomsobj.translate(dr)
            self.atoms.extend(atomsobj)
