# -*- coding: utf-8 -*-
"""
==============================================================================
Points class (:mod:`sknano.core.math._points`)
==============================================================================

.. currentmodule:: sknano.core.math._points

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import UserList
#from sknano.utils.geometric_shapes import Cuboid  # , Rectangle

__all__ = ['Points']


class Points(UserList):
    """Class for collection of `Point` objects.

    Parameters
    ----------
    points : {None, sequence, `Points`}, optional
        if not `None`, then a list of `Point` instance objects or an
        existing `Points` instance object.
    copylist : bool, optional
        perform shallow copy of points list
    deepcopy : bool, optional
        perform deepcopy of points list


    """

    def __init__(self, points=None, copylist=True, deepcopy=False):
        super(Points, self).__init__(initlist=points,
                                     copylist=copylist,
                                     deepcopy=deepcopy)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `Points`."""
        return "Points(points={!r})".format(self.data)

    def sort(self, key=None, reverse=False):

        if key is None:
            self.data.sort(key=attrgetter('x', 'y', 'z'),
                           reverse=reverse)
        else:
            self.data.sort(key=key, reverse=reverse)

    @property
    def x(self):
        """Return :math:`x` coordinates of `Point` objects as array."""
        return [point.x for point in self]

    @property
    def y(self):
        """Return :math:`y` coordinates of `Point` objects as array."""
        return [point.y for point in self]

    @property
    def z(self):
        """Return :math:`z` coordinates of `Point` objects as array."""
        return [point.z for point in self]

    def filter(self, condition, invert=False):
        """Filter `Points` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
            Boolean index array having same shape as the initial dimensions
            of the list of `Points` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

        Returns
        -------
        filtered_points : `Points`
            If `invert` is `False`, return the elements where `condition`
            is `True`.

            If `invert` is `True`, return the elements where `~condition`
            (i.e., numpy.invert(condition)) is `True`.

        """
        if invert:
            condition = ~condition
        return self.__class__(points=np.asarray(self)[condition].tolist())

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [point.rezero(epsilon=epsilon) for point in self]

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               deg2rad=False, transform_matrix=None):
        """Rotate point coordinates about arbitrary axis.

        Parameters
        ----------
        angle : float

        """
        [point.rotate(angle=angle, rot_axis=rot_axis,
                      anchor_point=anchor_point, deg2rad=deg2rad,
                      transform_matrix=transform_matrix) for point in self]

    def translate(self, t):
        """Translate point coordinates.

        Parameters
        ----------
        t : array_like
            3-elment array of :math:`x,y,z` components of translation vector
        """
        [point.translate(t) for point in self]
