# -*- coding: utf-8 -*-
"""
===============================================================================
misc core functions, preset data structures, etc. (:mod:`sknano.core._extras`)
===============================================================================

.. currentmodule:: sknano.core._extras

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['components', 'dimensions', 'xyz', 'xyz_axes', 'AttrDict']

components = dimensions = xyz = xyz_axes = ('x', 'y', 'z')


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
