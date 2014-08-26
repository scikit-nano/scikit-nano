# -*- coding: utf-8 -*-
"""
========================================================================
Helper functions (:mod:`sknano.utils.testing._funcs`)
========================================================================

.. currentmodule:: sknano.utils.testing._funcs

"""
from __future__ import absolute_import, division, print_function

#import numpy as np
from sknano.core.atoms import XAtom, XAtoms

__all__ = ['generate_atoms']


def generate_atoms(elements=None, structure=None, **kwargs):
    atoms = XAtoms()
    if elements is not None and isinstance(elements, list):
        for e in elements:
            atoms.append(XAtom(e))

    return atoms
