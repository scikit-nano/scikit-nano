# -*- coding: utf-8 -*-
"""
========================================================================
Helper functions (:mod:`sknano.utils.testing._funcs`)
========================================================================

.. currentmodule:: sknano.utils.testing._funcs

"""
from __future__ import absolute_import, division, print_function

#from inspect import getargspec
import importlib

#import numpy as np
from sknano.core.atoms import XAtom, XAtoms
from sknano.generators import STRUCTURE_GENERATORS

__all__ = ['generate_atoms']


def generate_atoms(elements=None, generator_class=None, **kwargs):
    if elements is not None and isinstance(elements, list):
        atoms = XAtoms()
        for e in elements:
            atoms.append(XAtom(e))
        return atoms
    elif generator_class is not None and \
            generator_class in STRUCTURE_GENERATORS:
        try:
            generator = getattr(importlib.import_module('sknano.generators'),
                                generator_class)
            structure = generator(**kwargs)
            return structure.atoms
        except ImportError as e:
            print(e)
