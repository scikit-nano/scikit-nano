# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for atoms. (:mod:`sknano.core.atoms.mixins`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins

Contents
========

Mixin `Atom`/`Atoms` classes
----------------------------

.. autosummary::
   :toctree: generated/

   PBCAtomsMixin

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._periodic_atoms import *

__all__ = [s for s in dir() if not s.startswith('_')]
