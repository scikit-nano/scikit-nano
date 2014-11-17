# -*- coding: utf-8 -*-
"""
===============================================================================
Fullerene generators (:mod:`sknano.generators._fullerene_generators`)
===============================================================================

.. currentmodule:: sknano.generators._fullerene_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from pkg_resources import resource_filename
import os
#import numpy as np

from sknano.io import XYZReader
from sknano.structures import Fullerene
#from sknano.utils.geometric_shapes import Cuboid
from ._base import GeneratorBase

__all__ = ['FullereneGenerator']


class FullereneGenerator(Fullerene, GeneratorBase):
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
    >>> fg.save_data(fname='C60.data')

    """

    def __init__(self, autogen=True, **kwargs):

        super(FullereneGenerator, self).__init__(**kwargs)

        if autogen:
            self.generate_structure_data()

    def generate_structure_data(self):
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
                #TODO: try to *intelligently* pick the best match
                CNfile = files[0]
            self.atoms = XYZReader(os.path.join(datadir, CNfile)).atoms

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True, **kwargs):
        """Save structure data.

        See :py:meth:`~sknano.generators.GeneratorBase.save_data` method
        for documentation.

        """
        if fname is None:
            fname = 'C{}'.format(self.N)

        super(FullereneGenerator, self).save_data(
            fname=fname, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM, **kwargs)
