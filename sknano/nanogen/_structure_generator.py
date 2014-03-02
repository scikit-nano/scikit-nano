# -*- coding: utf-8 -*-
"""
===================================================================
Structure generator (:mod:`sknano.nanogen._structure_generator`)
===================================================================

.. currentmodule:: sknano.nanogen._structure_generator

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext'

from ..chemistry import Atoms

__all__ = ['StructureGenerator']


class StructureGenerator(object):
    """Base class for structure generator classes."""

    def __init__(self):
        self._fname = None
        self._unit_cell = Atoms()
        self._structure_atoms = Atoms()

    @property
    def fname(self):
        """Structure file name."""
        return self._fname

    @property
    def unit_cell(self):
        """Return unit cell :py:class:`~sknano.chemistry.Atoms`."""
        return self._unit_cell

    @property
    def structure_atoms(self):
        """Return structure :py:class:`~sknano.chemistry.Atoms`."""
        return self._structure_atoms
