# -*- coding: utf-8 -*-
"""
========================================================================
Helper functions (:mod:`sknano.testing._funcs`)
========================================================================

.. currentmodule:: sknano.testing._funcs

"""
from __future__ import absolute_import, division, print_function

#from inspect import getargspec
import importlib

#import numpy as np
from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.core.refdata import element_symbols as periodic_table_of_elements
from sknano.generators import STRUCTURE_GENERATORS

__all__ = ['generate_atoms']


def generate_atoms(elements=None, generator_class=None, **kwargs):
    if elements is not None:
        if isinstance(elements, list):
            atoms = Atoms()
            for e in elements:
                atoms.append(Atom(element=e))
            return atoms
        elif isinstance(elements, str) and elements == 'periodic_table':
            atoms = Atoms()
            for e in periodic_table_of_elements:
                atoms.append(Atom(element=e))
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
