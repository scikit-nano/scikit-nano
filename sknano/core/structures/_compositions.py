# -*- coding: utf-8 -*-
"""
==============================================================================
Base composition classes (:mod:`sknano.core.structures._compositions`)
==============================================================================

.. currentmodule:: sknano.core.structures._compositions

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter

import numpy as np

from sknano.core import BaseClass, UserList
from sknano.core.math import convert_condition_str, rotation_matrix

__all__ = ['Composition', 'Compositions']


@total_ordering
class Composition(BaseClass):
    """Base class for abstract representation of a composition.

    Parameters
    ----------
    components : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instances or an
        `Atoms` instance.

    """
    def __init__(self, components=None, **kwargs):
        super().__init__(**kwargs)
        self.components = components
        self.fmtstr = "components={components!r}"

    def __eq__(self, other):
        """Test equality of two `Composition` object instances."""
        return self is other or self.components == other.components

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        return self.components < other.components

    def rezero(self, *args, **kwargs):
        assert not hasattr(super(), 'rezero')

    def rotate(self, **kwargs):
        assert not hasattr(super(), 'rotate')

    def translate(self, *args, **kwargs):
        assert not hasattr(super(), 'translate')

    def todict(self):
        """Return :class:`~python:dict` of `Composition` constructor \
            parameters."""
        return dict(components=self.components)


class Compositions(UserList):
    """Base class for collection of `Composition` objects.

    Parameters
    ----------
    compositions : {None, sequence, `Compositions`}, optional
        if not `None`, then a list of `Composition` instance objects or an
        existing `Compositions` instance object.

    """
    def __init__(self, compositions=None, update_item_class=True, **kwargs):
        super().__init__(initlist=compositions)

        self.fmtstr = "{compositions!r}"

    @property
    def __item_class__(self):
        return Composition

    def sort(self, key=attrgetter('id'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def Ncompositions(self):
        """Number of compositions in `Compositions`."""
        return len(self)

    @property
    def masses(self):
        """Return list of `Composition` masses."""
        return np.asarray([composition.mass for composition in self])

    def filter(self, condition, invert=False):
        """Filter `Compositions` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
        invert : bool, optional

        Returns
        -------
        filtered_compositions : `Compositions`

        """
        if isinstance(condition, str):
            condition = convert_condition_str(self, condition)
        if invert:
            condition = ~condition
        try:
            self.data = np.asarray(self)[condition].tolist()
        except AttributeError:
            self.data = np.asarray(self)[condition]

    def filtered(self, condition, invert=False):
        if isinstance(condition, str):
            condition = convert_condition_str(self, condition)
        if invert:
            condition = ~condition
        try:
            return self.__class__(components=np.asarray(self)[condition].tolist(),
                                  **self.kwargs)
        except AttributeError:
            return self.__class__(components=np.asarray(self)[condition],
                                  **self.kwargs)

    def get_compositions(self, asarray=False):
        """Return list of `Compositions`.

        Parameters
        ----------
        asarray : bool, optional

        Returns
        -------
        sequence or ndarray

        """
        if asarray:
            return np.asarray(self.data)
        else:
            return self.data

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [composition.rezero(epsilon=epsilon) for composition in self]

    def rotate(self, **kwargs):
        """Rotate `Composition` position vectors.

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
        if kwargs.get('transform_matrix', None) is None:
            kwargs['transform_matrix'] = rotation_matrix(**kwargs)
        [composition.rotate(**kwargs) for composition in self]

    def translate(self, t, fix_anchor_points=True):
        """Translate `Composition` position vectors by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [composition.translate(t, fix_anchor_point=fix_anchor_points)
         for composition in self]

    def todict(self):
        return dict(compositions=self.data)
