# -*- coding: utf-8 -*-
"""
===============================================================================
Fullerene generators (:mod:`sknano.generators._fullerene_generators`)
===============================================================================

.. currentmodule:: sknano.generators._fullerene_generators

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from pkg_resources import resource_filename
import os
# import numpy as np

from sknano.io import XYZReader
from sknano.core.structures import Fullerene
# from sknano.core.geometric_regions import Cuboid
from ._base import GeneratorBase

__all__ = ['FullereneGenerator']


class FullereneGenerator(GeneratorBase, Fullerene):
    """Fullerene structure generator class.

    Parameters
    ----------

    Raises
    ------

    Examples
    --------
    First, load the :class:`~sknano.generators.FullereneGenerator` class.

    >>> from sknano.generators import FullereneGenerator
    >>> fg = FullereneGenerator(60)
    >>> fg.save(fname='C60.data')

    """
    def generate(self):
        """Generate structure data."""
        self._atoms = XYZReader(self.datafile).atoms

    @classmethod
    def generate_fname(cls, datafile):
        return os.path.splitext(os.path.basename(datafile))[0]

    def save(self, fname=None, outpath=None, structure_format=None,
             center_centroid=True, **kwargs):
        """Save structure data.

        See :py:meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(self.datafile)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_centroid=center_centroid, **kwargs)
