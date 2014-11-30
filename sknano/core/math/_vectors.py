# -*- coding: utf-8 -*-
"""
==============================================================================
Vectors class (:mod:`sknano.core.math._vectors`)
==============================================================================

.. currentmodule:: sknano.core.math._vectors

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import UserList
from ._transforms import transformation_matrix
#from sknano.utils.geometric_shapes import Cuboid  # , Rectangle

__all__ = ['Vectors']


class Vectors(UserList):
    """Class for collection of `Vector` objects.

    Parameters
    ----------
    vectors : {None, sequence, `Vectors`}, optional
        if not `None`, then a list of `Vector` instance objects or an
        existing `Vectors` instance object.
    copylist : bool, optional
        perform shallow copy of vectors list
    deepcopy : bool, optional
        perform deepcopy of vectors list


    """

    def __init__(self, vectors=None, copylist=True, deepcopy=False):
        super(Vectors, self).__init__(initlist=vectors,
                                      copylist=copylist,
                                      deepcopy=deepcopy)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `Vectors`."""
        return "Vectors(vectors={!r})".format(self.data)

    def sort(self, key=None, reverse=False):

        if key is None:
            self.data.sort(key=attrgetter('p0', 'length'),
                           reverse=reverse)
        else:
            self.data.sort(key=key, reverse=reverse)

    @property
    def norms(self):
        """Return `Vector` :attr:`Vector.norm`\ s as array."""
        return [vector.norm for vector in self]

    def filter(self, condition, invert=False):
        """Filter `Vectors` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
            Boolean index array having same shape as the initial dimensions
            of the list of `Vectors` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

        Returns
        -------
        filtered_vectors : `Vectors`
            If `invert` is `False`, return the elements where `condition`
            is `True`.

            If `invert` is `True`, return the elements where `~condition`
            (i.e., numpy.invert(condition)) is `True`.

        """
        if invert:
            condition = ~condition
        return self.__class__(vectors=np.asarray(self)[condition].tolist())

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [vector.rezero(epsilon=epsilon) for vector in self]

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               deg2rad=False, transform_matrix=None, verbose=False):
        """Rotate `Vector`\ s coordinates.

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
        [vector.rotate(transform_matrix=transform_matrix) for vector in self]

    def translate(self, t, fix_anchor_points=False):
        """Translate `Vector`\ s by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [vector.translate(t, fix_anchor_point=fix_anchor_points)
         for vector in self]
