# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for MD analysis (:mod:`sknano.core.atoms._md_atom`)
===============================================================================

An `Atom` class for molecular dynamics structure analysis.

.. currentmodule:: sknano.core.atoms._md_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import numbers

from ._kdtree_atom import KDTAtom

__all__ = ['MDAtom']


class MDAtom(KDTAtom):
    """An `Atom` class for molecular dynamics structure analysis.

    Parameters
    ----------

    """
    def __init__(self, **kwargs):
        super(MDAtom, self).__init__(**kwargs)

    def __str__(self):
        """Return a nice string representation of `MDAtom`."""
        return super(MDAtom, self).__str__()

    def __repr__(self):
        """Return canonical string representation of `MDAtom`."""
        #strrep = "Atom(element={element!r}, atomID={atomID!r}, " + \
        #    "moleculeID={moleculeID!r}, atomtype={atomtype!r}, " + \
        #    "q={q!r}, m={m!r}, x={x:.6f}, y={y:.6f}, z={z:.6f}, " + \
        #    "CN={CN!r}, NN={NN!r})"
        #parameters = dict(element=self.element, atomID=self.atomID,
        #                  moleculeID=self.moleculeID, atomtype=self.atomtype,
        #                  q=self.q, m=self.m, x=self.x, y=self.y, z=self.z,
        #                  CN=self.CN, NN=self.NN)
        #return strrep.format(**parameters)
        return super(MDAtom, self).__repr__()
