# -*- coding: utf-8 -*-
"""
========================================================================
Class for XYZ atom (:mod:`sknano.structure_io.atoms._xyz_atom`)
========================================================================

.. currentmodule:: sknano.structure_io.atoms._xyz_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._atom import Atom

__all__ = ['XYZAtom']


class XYZAtom(Atom):
    """`Atom` class for `xyz` format.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    x, y, z : float, optional
        :math:`x, y, z` components of `Atom` position vector relative to
        origin.

    """

    def __init__(self, element=None, x=None, y=None, z=None):

        super(XYZAtom, self).__init__(element=element, x=x, y=y, z=z)
