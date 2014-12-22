# -*- coding: utf-8 -*-
"""
===============================================================================
MWNT structure generator (:mod:`sknano.generators._mwnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators._mwnt_generator

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

#import copy

import numpy as np

#from sknano.core import pluralize
#from sknano.core.math import Vector
from sknano.structures import MWNT
from sknano.utils.geometric_shapes import Cuboid
from ._base import GeneratorBase
from ._swnt_generator import SWNTGenerator

__all__ = ['MWNTGenerator']


class MWNTGenerator(MWNT, GeneratorBase):
    """Class for generating single, `MWNT`.

    .. versionchanged:: 0.2.20

       `MWNTGenerator` no longer generates MWNT *bundles*, only *single*
       MWNTs. To generate bundled MWNT structure data, use the
       `MWNTBundleGenerator` class.

    .. versionadded:: 0.2.8

    Parameters
    ----------
    Ch_list : :class:`python:list`, optional
        (:attr:`~SWNT.n`, :attr:`~SWNT.m`) for each `SWNT` wall in `MWNT`.
    Nwalls : int, optional
        Number of `SWNT` walls in `MWNT`.
    Lz : float, optional
        `MWNT` length in **nanometers**.
    min_shell_diameter : float, optional
        Minimum `MWNT` wall diameter, in units of **Angstroms**.
    max_shell_diameter : float, optional
        Maximum `MWNT` wall diameter, in units of **Angstroms**.
    max_shells : int, optional
        Maximum number of `MWNT` walls.
    chiral_types : {None, 'armchair', 'zigzag', 'achiral', 'chiral'}, optional
        If `None`, the :attr:`~SWNT.chiral_type` of each `MWNT` walls
        will be random and determined by the set of randomly selected
        chiral indices (:attr:`~SWNT.n`, :attr:`~SWNT.m`).
    shell_spacing : float, optional
        Inter-wall spacing in units of **Angstroms**.
        Default value is the van der Waals interaction distance of 3.4
        Angstroms.
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~MWNTGenerator.generate_structure_data`.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom` 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms, in units of **Angstroms**.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.generators import MWNTGenerator
    >>> mwnt = MWNTGenerator(Nwalls=5, max_shell_diameter=50,
    ...                      Lz=1.0, fix_Lz=True)
    >>> mwnt.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_1cellx1cellx4.06cells-01.png

    """
    def __init__(self, autogen=True, **kwargs):

        super(MWNTGenerator, self).__init__(**kwargs)

        if autogen:
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data.

        .. todo::

           Load the diameter and chirality data from file instead of
           generating it every time.

        """

        self.structure_data.clear()
        for swnt in self.shells:
            self.atoms.extend(SWNTGenerator(**swnt.todict()).atoms)

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  anchor_point=None, center_CM=True, savecopy=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save_data` method
        for documentation.

        """
        if fname is None:
            Nshells = '{}shell_mwnt'.format(self.Nshells)
            chiralities = '@'.join([str(Ch).replace(' ', '') for
                                    Ch in self.Ch_list])

            fname = '_'.join((Nshells, chiralities))

        Lz_cutoff = (10 * self.L0 + 1) / 2
        pmax = np.asarray([np.inf, np.inf, Lz_cutoff])
        pmin = -pmax
        region_bounds = Cuboid(pmin=pmin.tolist(), pmax=pmax.tolist())
        region_bounds.update_region_limits()

        self.atoms.clip_bounds(region_bounds, center_before_clipping=True)

        if center_CM:
            self.atoms.center_CM()

        super(MWNTGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=False,
            savecopy=savecopy, **kwargs)
