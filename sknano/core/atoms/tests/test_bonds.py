#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename
import numpy as np
#from sknano.core.atoms import Bond, Bonds
from sknano.generators import SWNTGenerator
from sknano.io import DATAReader


def test1():
    infile = resource_filename('sknano', 'data/nanotubes/1005_5cells.data')
    data = DATAReader(infile)
    atoms = data.atoms
    #atoms._update_nearest_neighbors()
    #atoms.update_pyramidalization_angles()
    bonds = atoms.bonds
    assert_equal(bonds.Nbonds, atoms.coordination_numbers.sum())
    print(bonds.mean_length)
    print(bonds.atoms.Natoms)
    #mean_angle = bonds.mean_angle

    #print('bonds.mean_angle: {}'.format(bonds.mean_angle))


def test2():
    atoms = SWNTGenerator(n=10, m=10, nz=3).atoms
    #assert_equal(atoms.Natoms, 80)
    atoms.assign_unique_ids()
    #print(atoms.bonds.mean_length)
    #print(atoms.bonds.mean_angle)
    misalignment_angles = atoms.poav_misalignment_angles
    print(len(misalignment_angles))
    print(np.degrees(misalignment_angles))
    print(np.degrees(atoms.pyramidalization_angles))


def test3():
    atoms = SWNTGenerator(n=20, m=0, nz=2).atoms
    assert_equal(atoms.Natoms, 160)
    atoms.assign_unique_ids()
    atoms._update_nearest_neighbors()
    #print(atoms.bonds.mean_length)
    #print(atoms.bonds.mean_angle)
    atom0 = atoms[0]
    atom0bonds = atom0.bonds
    print('atom0: {}'.format(atom0))
    print('atom0.r: {}'.format(atom0.r))
    for NN in atom0.NN:
        print('NN.r: {}'.format(NN.r))
    print(atom0bonds.atoms.Natoms)
    print(atom0bonds.atoms.CM)
    print('atoms.bonds.Nbonds: {}'.format(atoms.bonds.Nbonds))
    print('atoms.bonds.atoms.Natoms: {}'.format(atoms.bonds.atoms.Natoms))


if __name__ == '__main__':
    nose.runmodule()
