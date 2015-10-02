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

from pkg_resources import resource_filename
import os
# import numpy as np

from sknano.io import XYZReader
from sknano.structures import Fullerene
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
    >>> fg = FullereneGenerator(N=60)
    >>> fg.save(fname='C60.data')

    """
    def generate(self):
        """Generate structure data."""
        CNdir = 'C' + str(self.N)
        CNfile = 'C' + str(self.N)
        if self.PG is not None:
            CNfile += '-' + self.PG
        if self.Ni is not None:
            CNfile += '-{}'.format(self.Ni)
        CNfile += '.xyz'

        datadir = os.path.join(resource_filename('sknano', 'data/fullerenes'),
                               CNdir)
        files = os.listdir(datadir)
        if len(files) > 0:
            if CNfile not in files:
                # TODO: try to *intelligently* pick the best match
                CNfile = files[0]
            self._atoms = XYZReader(os.path.join(datadir, CNfile)).atoms

    @classmethod
    def generate_fname(cls, N):
        return 'C{}'.format(N)

    def save(self, fname=None, outpath=None, structure_format=None,
             center_centroid=True, **kwargs):
        """Save structure data.

        See :py:meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(self.N)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_centroid=center_centroid, **kwargs)
