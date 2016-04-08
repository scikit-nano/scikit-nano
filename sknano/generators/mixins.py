# -*- coding: utf-8 -*-
"""
==============================================================================
Generator mixin classes (:mod:`sknano.generators.mixins`)
==============================================================================

.. currentmodule:: sknano.generators.mixins

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import copy

__all__ = ['CappedNanotubeGeneratorMixin']


class CappedNanotubeGeneratorMixin:
    """Mixin class for generating capped nanotubes."""

    def generate_endcaps(self):
        pass
