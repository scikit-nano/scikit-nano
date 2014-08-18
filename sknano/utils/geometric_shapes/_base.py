# -*- coding: utf-8 -*-
"""
===============================================================================
Base classes for geometric regions (:mod:`sknano.utils.geometric_shapes._base`)
===============================================================================

.. currentmodule:: sknano.utils.geometric_shapes._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod, abstractproperty

__all__ = ['GeometricRegion']


class GeometricRegion(object):
    """Abstract base class for geometric regions."""
    __metaclass__ = ABCMeta

    @abstractproperty
    def centroid(self):
        """Centroid of 2D geometric region."""
        raise NotImplementedError

    @abstractmethod
    def contains_point(self):
        """Check if point is contained within geometric region."""
        raise NotImplementedError
