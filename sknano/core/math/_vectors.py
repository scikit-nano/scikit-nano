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

    def rotate(self, angle, rot_axis=None, anchor_vector=None,
               deg2rad=False):
        """Rotate vector coordinates about arbitrary axis.

        Parameters
        ----------
        angle : float

        """
        [vector.rotate(angle, rot_axis=rot_axis, anchor_vector=anchor_vector,
                       deg2rad=deg2rad) for vector in self]

    def translate(self, t):
        """Translate vector coordinates.

        Parameters
        ----------
        t : array_like
            3-elment array of :math:`x,y,z` components of translation vector
        """
        [vector.translate(t) for vector in self]
