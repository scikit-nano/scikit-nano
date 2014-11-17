# -*- coding: utf-8 -*-
"""
==============================================================================
Structure analysis module (:mod:`sknano.utils.analysis._structure_analyzer`)
==============================================================================

.. currentmodule:: sknano.utils.analysis._structure_analyzer

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['StructureAnalyzer']

import importlib

import numpy as np

from sknano.core.math import vector as vec


class StructureAnalyzer(object):
    u"""Class for structure analysis.

    Parameters
    ----------
    structure_data : str
        structure data file

    """
    def __init__(self, structure_data, **kwargs):
        self.structure_data = structure_data
        self.structure_data.atoms.update_attrs()

    def analyze_POAV(self):
        for atom in self.structure_data.atoms:
            # the central atom must have 3 bonds
            if atom.bonds.Nbonds == 3:
                for POAV_name in ('POAV1', 'POAV2', 'POAVR'):
                    POAV_class = \
                        getattr(
                            importlib.import_module('sknano.utils.analysis'),
                            POAV_name)
                    setattr(atom, POAV_name, POAV_class(atom.bonds))
                    POAV = getattr(atom, POAV_name)
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
                        if getattr(NN, POAV_name) is not None:
                            NN_POAV = getattr(NN, POAV_name)
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
                            misalignment_angles.append(np.abs(
                                np.pi / 2 - vec.angle(NN_POAV.Vpi, nvec)))
                        else:
                            misalignment_angles.append(np.nan)

                    POAV.pyramidalization_angles = pyramidalization_angles
                    POAV.misalignment_angles = misalignment_angles
                    POAV.sigma_pi_angles = sigma_pi_angles
