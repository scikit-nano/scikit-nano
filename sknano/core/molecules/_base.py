# -*- coding: utf-8 -*-
"""
==============================================================================
Base classes for molecules package (:mod:`sknano.core.molecules._base`)
==============================================================================

.. currentmodule:: sknano.core.molecules._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from sknano.core import UserList

__all__ = ['MoleculeList']


class MoleculeList(UserList):
    """Base class for collection of `Molecule` objects.

    Parameters
    ----------
    molecules : {None, sequence, `molecules`}, optional
        if not `None`, then a list of `Molecule` instance objects or an
        existing `molecules` instance object.
    copylist : bool, optional
        perform shallow copy of molecules list
    deepcopy : bool, optional
        perform deepcopy of molecules list

    """

    def __init__(self, molecules=None, copylist=True, deepcopy=False):
        super(MoleculeList, self).__init__(initlist=molecules,
                                           copylist=copylist,
                                           deepcopy=deepcopy)
