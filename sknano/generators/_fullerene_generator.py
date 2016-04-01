# -*- coding: utf-8 -*-
"""
===============================================================================
Fullerene generator class (:mod:`sknano.generators._fullerene_generator`)
===============================================================================

.. currentmodule:: sknano.generators._fullerene_generator

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from pkg_resources import resource_filename
import os
# import numpy as np

from sknano.core.structures import Fullerene
# from sknano.core.geometric_regions import Cuboid
from ._base import NanoStructureGenerator, GeneratorMixin

__all__ = ['FullereneGenerator']


class FullereneGenerator(GeneratorMixin, NanoStructureGenerator, Fullerene):
    """Fullerene structure generator class.

    Parameters
    ----------

    Raises
    ------

    Examples
    --------
    First, load the :class:`~sknano.generators.FullereneGenerator` class.

    >>> from sknano.generators import FullereneGenerator
    >>> buckyball = FullereneGenerator(60)
    >>> buckyball.save()

    .. image:: /images/buckyball-1.png

    """

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
