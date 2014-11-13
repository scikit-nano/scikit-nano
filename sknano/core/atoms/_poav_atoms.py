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

import numpy as np

from sknano.core.math import vector as vec
from ._kdtree_atoms import KDTAtoms
from ._poav_atom import POAV1, POAV2, POAVR

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
    _atomattrs = KDTAtoms._atomattrs + ['POAVlist', 'POAV1', 'POAV2', 'POAVR']

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        super(POAVAtoms, self).__init__(atoms=atoms,
                                        copylist=copylist,
                                        deepcopy=deepcopy)

    def update_attrs(self):
        super(POAVAtoms, self).update_attrs()
        self._update_POAVlist()

    def _update_POAVlist(self):
        """Compute `POAV1`, `POAV2`, `POAVR`."""
        for atom in self:
            # the central atom must have 3 bonds
            if atom.bonds.Nbonds == 3:
                atom.POAVlist = [POAV1(atom.bonds), POAV2(atom.bonds),
                                 POAVR(atom.bonds)]
                for i, POAV in enumerate(atom.POAVlist):
                    pyramidalization_angles = []
                    misalignment_angles = []
                    sigma_pi_angles = []
                    for bond, NN in zip(atom.bonds, atom.NN):
                        # first compute the pyramidalization angle
                        sigma_pi_angle = vec.angle(POAV.Vpi, bond.vector)
                        if sigma_pi_angle < np.pi / 2:
                            sigma_pi_angle = np.pi - sigma_pi_angle
                        sigma_pi_angles.append(sigma_pi_angle)
                        pyramidalization_angles.append(
                            sigma_pi_angle - np.pi / 2)

                        # the bonded atom must have a POAV to compute the
                        # misalignment angles
                        if NN.POAVlist is not None:
                            # compute vector that is orthogonal to the plane
                            # defined by the bond vector and the POAV of the
                            # center atom.
                            nvec = vec.cross(bond.vector, POAV.Vpi)

                            # the misalignment angle is the angle between the
                            # nearest neighbor's POAV and the plane defined by
                            # the bond vector and the POAV of the center atom,
                            # which is pi/2 minus the angle between
                            # the NN POAV and the normal vector to the plane
                            # computed above.
                            misalignment_angles.append(
                                np.abs(np.pi / 2 -
                                       vec.angle(NN.POAVlist[i].Vpi, nvec)))
                        else:
                            misalignment_angles.append(np.nan)

                    POAV.pyramidalization_angles = pyramidalization_angles
                    POAV.misalignment_angles = misalignment_angles
                    POAV.sigma_pi_angles = sigma_pi_angles
