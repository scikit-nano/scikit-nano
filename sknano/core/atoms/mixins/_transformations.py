# -*- coding: utf-8 -*-
"""
=====================================================================================
Mixin classes for transformations (:mod:`sknano.core.atoms.mixins._transformations`)
=====================================================================================

.. currentmodule:: sknano.core.atoms.mixins._transformations

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import numpy as np
# import pandas as pd

from sknano.core.math import transformation_matrix

__all__ = ['AtomTransformationsMixin', 'AtomsTransformationsMixin']


class AtomTransformationsMixin:
    """Mixin `Atom` class for performing affine transformations."""

    def rotate(self, **kwargs):
        """Rotate `Atom` position vector.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        transform_matrix = kwargs.get('transform_matrix', None)
        if transform_matrix is None:
            kwargs['transform_matrix'] = transformation_matrix(**kwargs)

        try:
            self.lattice.rotate(**kwargs)
        except AttributeError:
            pass
        self.r.rotate(**kwargs)
        self.r0.rotate(**kwargs)

    def translate(self, t, fix_anchor_point=True, cartesian=True):
        """Translate :class:`Atom` position vector by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : :class:`~python:bool`, optional
        cartesian : :class:`~python:bool`, optional

        """
        try:
            if not cartesian:
                t = self.lattice.fractional_to_cartesian(t)
            # if not fix_anchor_point:
            #     self.lattice.translate(t)
            self.lattice.translate(t)
        except AttributeError:
            pass

        # TODO compare timing benchmarks for translation of position vector
        self.r.translate(t, fix_anchor_point=fix_anchor_point)
        self.r0.translate(t, fix_anchor_point=fix_anchor_point)
        # self.r += t


class AtomsTransformationsMixin:
    """Mixin `Atoms` class for performing affine transformations."""

    def rotate(self, **kwargs):
        """Rotate `Atom` vectors.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        if kwargs.get('anchor_point', None) is None:
            try:
                kwargs['anchor_point'] = self.lattice.offset
            except AttributeError:
                try:
                    kwargs['anchor_point'] = self.bounding_region.centroid
                except AttributeError:
                    kwargs['anchor_point'] = self.centroid

        transform_matrix = kwargs.get('transform_matrix', None)
        if transform_matrix is None:
            kwargs['transform_matrix'] = transformation_matrix(**kwargs)
        [atom.rotate(**kwargs) for atom in self]

    def translate(self, t, fix_anchor_points=True, cartesian=True):
        """Translate `Atom` vectors by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [atom.translate(t, fix_anchor_point=fix_anchor_points,
                        cartesian=cartesian) for atom in self]
