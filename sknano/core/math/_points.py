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
from ._transforms import transformation_matrix
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
               rot_point=None, from_vector=None, to_vector=None,
               deg2rad=False, transform_matrix=None, verbose=False):
        """Rotate `Point`\ s coordinates.

        Parameters
        ----------
        angle : float
        rot_axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        deg2rad : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, rot_axis=rot_axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, deg2rad=deg2rad,
                                      verbose=verbose)
        [point.rotate(transform_matrix=transform_matrix) for point in self]

    def translate(self, t):
        """Translate `Point`\ s by :class:`Vector` `t`.

        Parameters
        ----------
        v : :class:`Vector`

        """
        [point.translate(t) for point in self]
