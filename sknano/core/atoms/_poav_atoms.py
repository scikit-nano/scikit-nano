# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for POAV analysis (:mod:`sknano.core.atoms._poav_atoms`)
===============================================================================

An `Atoms` class for POAV analysis

.. currentmodule:: sknano.core.atoms._poav_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import numbers
#import warnings
#warnings.filterwarnings('ignore', "Mean of empty slice.")
#warnings.filterwarnings('ignore',
#                        'invalid value encountered in double_scalars')

import numpy as np

from sknano.core.math import Vector, vector as vec
#from ._bond import Bond
#from ._bonds import Bonds
from ._kdtree_atoms import KDTAtoms

__all__ = ['POAVAtoms']


class POAVAtoms(KDTAtoms):
    """An `Atoms` sub-class for POAV analysis.

    Sub-class of `KDTAtoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.POAVAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `POAVAtoms`}, optional
        if not `None`, then a list of `POAVAtom` instance objects or an
        existing `POAVAtoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list

    """
    _atomattrs = KDTAtoms._atomattrs + \
        ['pyramidalization_angle', 'sigma_bond_angle', 'poav', 'poma']

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        super(POAVAtoms, self).__init__(atoms=atoms,
                                        copylist=copylist,
                                        deepcopy=deepcopy)

    #@property
    #def mean_nonzero_poma(self):
    #    #self._update_poma()
    #    return np.ma.mean(np.ma.fix_invalid([np.ma.mean(np.ma.compressed(
    #        np.ma.masked_values(np.ma.fix_invalid(atom.poma), 0)))
    #        for atom in self]))

    #@property
    #def nonzero_poma(self):
    #    #self._update_poma()
    #    return [np.ma.masked_values(np.ma.fix_invalid(atom.poma), 0)
    #            for atom in self]

    @property
    def poma(self):
        """Return per-atom list of POAV misalignment angles."""
        #self._update_poma()
        return np.ma.asarray([np.ma.fix_invalid(atom.poma) for atom in self])

    @property
    def pyramidalization_angles(self):
        #self._update_pyramidalization_angles()
        angles = []
        [angles.append(atom.pyramidalization_angle) for atom in self if
         atom.poav is not None]
        return np.asarray(angles)

    def update_attrs(self):
        super(POAVAtoms, self).update_attrs()
        self._update_pyramidalization_angles()
        self._update_poma()

    def _update_poma(self):
        #self._update_pyramidalization_angles()
        for atom in self:
            poma = []
            for i, NN in enumerate(atom.NN):
                bond = atom.bonds[i]
                if atom.poav is not None and NN.poav is not None:
                    nvec = vec.cross(bond.vector, atom.poav)
                    poma.append(np.abs(np.pi / 2 - vec.angle(NN.poav, nvec)))
                else:
                    poma.append(np.nan)
            atom.poma = poma

    def _update_pyramidalization_angles(self):
        #self._update_bonds()
        for atom in self:
            if atom.bonds.Nbonds == 3:
                b1, b2, b3 = atom.bonds
                v21 = Vector(b2.vector - b1.vector, p0=b1.vector.p)
                v31 = Vector(b3.vector - b1.vector, p0=b1.vector.p)
                poav = vec.cross(v21, v31)
                atom.poav = poav.unit_vector
                atom.sigma_bond_angle = vec.angle(atom.poav, b1.vector)
                if atom.sigma_bond_angle < np.pi / 2:
                    atom.sigma_bond_angle = np.pi - atom.sigma_bond_angle
                    atom.poav = -atom.poav
                atom.pyramidalization_angle = atom.sigma_bond_angle - np.pi / 2
