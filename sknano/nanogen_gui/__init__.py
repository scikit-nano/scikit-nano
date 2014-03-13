# -*- coding: utf-8 -*-
"""
===================================================
NanoGen GUI front-ends (:mod:`sknano.nanogen_gui`)
===================================================

.. currentmodule:: sknano.nanogen_gui

.. versionadded:: 0.2.24
   `NanoGen` MVC GUI implemented.

.. seealso:: CLI module :py:mod:`sknano.scripts.nanogen`

Contents
========

.. autosummary::
   :toctree: generated/

   NanoGen

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext en'

from .nanogenui import NanoGen

__all__ = [s for s in dir() if not s.startswith('_')]
