# -*- coding: utf-8 -*-
"""
==============================================================================
Base classes for atoms package (:mod:`sknano.core.atoms._base`)
==============================================================================

.. currentmodule:: sknano.core.atoms._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from sknano.core import UserList

__all__ = ['AtomList', 'BondList']


class AtomList(UserList):
    """Base class for collection of `Atom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects or an
        existing `Atoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        super(AtomList, self).__init__(initlist=atoms,
                                       copylist=copylist,
                                       deepcopy=deepcopy)


class BondList(UserList):
    """Base class for collection of atom `Bonds`.

    Parameters
    ----------
    bonds : {None, sequence, `Bonds`}, optional
        if not `None`, then a list of `Bond` instance objects
    copylist : bool, optional
        perform shallow copy of bonds list
    deepcopy : bool, optional
        perform deepcopy of bonds list

    """
    def __init__(self, bonds=None, copylist=True, deepcopy=False):
        super(BondList, self).__init__(initlist=bonds,
                                       copylist=copylist,
                                       deepcopy=deepcopy)
