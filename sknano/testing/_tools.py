# -*- coding: utf-8 -*-
"""
========================================================================
Test fixtures (:mod:`sknano.testing._tools`)
========================================================================

.. currentmodule:: sknano.testing._tools

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import os
import unittest

from pkg_resources import resource_filename

from sknano.core.geometric_regions import Parallelogram, Rectangle, Square, \
    Ellipse, Circle, Triangle, Parallelepiped, Cuboid, Cube, Ellipsoid, \
    Sphere, Cylinder, Cone

from sknano.io import DATAReader, DUMPReader, PDBReader, XYZReader, \
    DATAData, DUMPData, PDBData, XYZData

from ._funcs import generate_atoms, generate_structure

__all__ = ['AtomsTestFixture', 'TempfileTestFixture', 'GeneratorTestFixture',
           'IOTestFixture', 'DUMPTestFixture', 'GeometricRegionsTestFixture',
           'Geometric2DRegionsTestFixture', 'Geometric3DRegionsTestFixture']


class AtomsTestFixture(unittest.TestCase):
    """Mixin :class:`~python:unittest.TestCase` class.

    Defines properties which return structures or atoms objects.

    """
    @property
    def atoms(self):
        """:class:`~sknano.generators.SWNTGenerator` \
            :class:`~sknano.core.atoms.StructureAtoms`."""
        return self.swnt.atoms

    @property
    def structure(self):
        """:class:`~sknano.generators.SWNTGenerator`"""
        return self.swnt

    @property
    def buckyball(self):
        """:class:`~sknano.generators.FullereneGenerator`."""
        buckyball = generate_structure(generator_class='FullereneGenerator',
                                       N=60)
        buckyball.update_attrs()
        return buckyball

    @property
    def graphene(self):
        """:class:`~sknano.generators.GrapheneGenerator` structure."""
        graphene = generate_structure(generator_class='GrapheneGenerator',
                                      armchair_edge_length=10,
                                      zigzag_edge_length=10)
        graphene.update_attrs()
        return graphene

    @property
    def bilayer_graphene(self):
        """:class:`~sknano.generators.GrapheneGenerator` structure."""
        blg = generate_structure(generator_class='BilayerGrapheneGenerator',
                                 armchair_edge_length=10,
                                 zigzag_edge_length=10)
        blg.update_attrs()
        return blg

    @property
    def periodic_table(self):
        """:class:`~sknano.core.atoms.StructureAtoms`"""
        periodic_table = generate_atoms(elements='periodic_table')
        periodic_table.assign_unique_ids()
        periodic_table.assign_unique_types()
        return periodic_table

    @property
    def dumpdata1(self):
        """:class:`~sknano.io.DUMPReader` object."""
        dumpfile = \
            resource_filename('sknano',
                              'data/lammpstrj/' +
                              'irradiated_graphene.system.dump.02000')

        return DUMPReader(dumpfile,
                          dumpattrmap={'c_atom_pe': 'pe', 'c_atom_ke': 'ke'},
                          atomattrmap={('type', 'element'): {1: 'C', 2: 'Ar'}})

    @property
    def swnt(self):
        """:class:`~sknano.generators.SWNTGenerator` structure."""
        swnt = generate_structure(generator_class='SWNTGenerator',
                                  n=5, m=0, nz=5)
        swnt.update_attrs()
        return swnt


class TempfileTestFixture(unittest.TestCase):
    """Mixin :class:`~python:unittest.TestCase` class temp file I/O.

    Defines setUp/tearDown methods to keep track of and delete temporary files
    created by unit tests.

    Attributes
    ----------
    tmpdata : :class:`~python:list`

    """
    def setUp(self):
        """Initialize list for names of temporary files."""
        self.tmpdata = []

    def tearDown(self):
        """Remove files in :attr:`~TempfileTestFixture.tmpdata`."""
        for f in self.tmpdata:
            try:
                os.remove(f)
            except IOError:
                continue


class GeneratorTestFixture(TempfileTestFixture):
    """Mixin :class:`~python:unittest.TestCase` class for \
        :mod:`sknano.generators` tests."""
    pass


class IOTestFixture(TempfileTestFixture):
    """Mixin :class:`~python:unittest.TestCase` class for \
        :mod:`sknano.io` unit tests.

    Defines lazy properties to get IO test data and IO class instances.

    """
    @property
    def data_reader(self):
        """:class:`~sknano.io.DATAReader` instance."""
        datafile = \
            resource_filename('sknano', 'data/nanotubes/1010_1cell.data')
        return DATAReader(datafile)

    @property
    def dump_reader(self):
        """:class:`~sknano.io.DUMPReader` instance."""
        dumpfile = \
            resource_filename('sknano', 'data/lammpstrj/0500_29cells.dump')

        return DUMPReader(dumpfile,
                          dumpattrmap={'c_peratom_pe': 'pe',
                                       'c_peratom_ke': 'ke'},
                          atomattrmap={('type', 'element'): {1: 'C'}})

    @property
    def pdb_reader(self):
        """:class:`~sknano.io.PDBReader` instance."""
        pdbfile = resource_filename('sknano', 'data/nanotubes/1010_1cell.pdb')
        return PDBReader(pdbfile)

    @property
    def xyz_reader(self):
        """:class:`~sknano.io.XYZReader` instance."""
        xyzfile = resource_filename('sknano', 'data/nanotubes/1010_1cell.xyz')
        return XYZReader(xyzfile)

    @property
    def datadata(self):
        """:class:`~sknano.io.DATAData` instance."""
        datafile = \
            resource_filename('sknano', 'data/nanotubes/1010_1cell.data')
        return DATAData(datafile)

    @property
    def dumpdata(self):
        """:class:`~sknano.io.DUMPData` instance."""
        dumpfile = \
            resource_filename('sknano', 'data/lammpstrj/0500_29cells.dump')

        return DUMPData(dumpfile, dumpattrmap={'c_peratom_pe': 'pe',
                                               'c_peratom_ke': 'ke',
                                               'v_speed': 'v.magnitude'},
                        atomattrmap={('type', 'element'): {1: 'C'}})

    @property
    def pdbdata(self):
        """:class:`~sknano.io.PDBData` instance."""
        pdbfile = resource_filename('sknano', 'data/nanotubes/1010_1cell.pdb')
        return PDBData(pdbfile)

    @property
    def xyzdata(self):
        """:class:`~sknano.io.XYZData` instance."""
        xyzfile = resource_filename('sknano', 'data/nanotubes/1010_1cell.xyz')
        return XYZData(xyzfile)


class DUMPTestFixture(IOTestFixture):
    """Mixin :class:`~python:unittest.TestCase` class for \
        :class:`~sknano.io.DUMPIO` unit tests.
    """
    def print_dumpattrs(self, dump):
        """Print dump attributes."""
        print('DUMP:\n{}'.format(dump))


class GeometricRegionsTestFixture(unittest.TestCase):
    """Mixin :class:`~python:unittest.TestCase` class for \
        :mod:`sknano.core.geometric_regions` unit tests.
    """

    def setUp(self):
        """Initialize list of regions."""
        self.regions = []
        print()

    def print_regions(self):
        """Pretty print regions."""
        for r in self.regions:
            print('{}\n'.format(r))

    def tearDown(self):
        """Pretty print regions."""
        self.print_regions()
        print()


class Geometric2DRegionsTestFixture(GeometricRegionsTestFixture):
    """Mixin :class:`~python:unittest.TestCase` class for \
        :class:`sknano.core.geometric_regions.Geometric2DRegions` unit tests.
    """
    @property
    def parallelogram(self):
        """Return :class:`~sknano.core.geometric_regions.Parallelogram`"""
        return Parallelogram()

    @property
    def rectangle(self):
        """Return :class:`~sknano.core.geometric_regions.Rectangle`"""
        return Rectangle()

    @property
    def square(self):
        """Return :class:`~sknano.core.geometric_regions.Square`"""
        return Square()

    @property
    def ellipse(self):
        """Return :class:`~sknano.core.geometric_regions.Ellipse`"""
        return Ellipse()

    @property
    def circle(self):
        """Return :class:`~sknano.core.geometric_regions.Circle`"""
        return Circle()

    @property
    def triangle(self):
        """Return :class:`~sknano.core.geometric_regions.Triangle`"""
        return Triangle()


class Geometric3DRegionsTestFixture(GeometricRegionsTestFixture):
    """Mixin :class:`~python:unittest.TestCase` class for \
        :class:`sknano.core.geometric_regions.Geometric3DRegions` unit tests.
    """

    @property
    def parallelepiped(self):
        """Return :class:`~sknano.core.geometric_regions.Parallelepiped`"""
        return Parallelepiped()

    @property
    def cuboid(self):
        """Return :class:`~sknano.core.geometric_regions.Cuboid`"""
        return Cuboid()

    @property
    def cube(self):
        """Return :class:`~sknano.core.geometric_regions.Cube`"""
        return Cube()

    @property
    def ellipsoid(self):
        """Return :class:`~sknano.core.geometric_regions.Ellipsoid`"""
        return Ellipsoid()

    @property
    def sphere(self):
        """Return :class:`~sknano.core.geometric_regions.Sphere`"""
        return Sphere()

    @property
    def cylinder(self):
        """Return :class:`~sknano.core.geometric_regions.Cylinder`"""
        return Cylinder()

    @property
    def cone(self):
        """Return :class:`~sknano.core.geometric_regions.Cone`"""
        return Cone()
