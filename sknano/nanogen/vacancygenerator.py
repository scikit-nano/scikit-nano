# -*- coding: utf-8 -*-
"""
=============================================================================
Generate vacancies in structures (:mod:`sknano.nanogen.vacancygenerator`)
=============================================================================

Module for generating vacancies in nano-structures.

.. currentmodule:: sknano.nanogen.vacancygenerator

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

import itertools
import os

import numpy as np

from pkshared.tools.refdata import CCbond

from .graphene import GrapheneGenerator, GrapheneGeneratorError
from .nanotube import NanotubeBundleGenerator, NanotubeGeneratorError
from ..structure_io import DATAReader, DATAWriter, XYZWriter, \
    XYZ2DATAConverter, StructureFormatError, supported_structure_formats

__all__ = ['GrapheneVacancyGenerator',
           'NanotubeVacancyGenerator',
           'VacancyGenerator', 'VacancyGeneratorError']


class VacancyGeneratorError(Exception):
    """Base class for VacancyGenerator exceptions."""
    pass


class VacancyGenerator(object):
    """Generate vacancies in structure data.

    Parameters
    ----------
    fname : str
        structure data filename
    structure_format : {None, str}, optional
        chemical file format of saved structure data.
        If ``None``, then guess based on ``fname`` file extension.
        Otherwise, must be one of:

            - xyz
            - data

    verbose : bool, optional
        Verbose output

    """
    def __init__(self, fname=str, structure_format=None, verbose=False):

        self.verbose = verbose

        if fname.endswith(supported_structure_formats) and \
                structure_format is None:
            for ext in supported_structure_formats:
                if fname.endswith(ext):
                    structure_format = ext
                    break
        else:
            if (not fname.endswith(supported_structure_formats) and
                    structure_format is None) or \
                    (structure_format is not None and
                        structure_format not in supported_structure_formats):
                raise StructureFormatError(
                    '{} is not a supported structure format'.format(
                        structure_format))

        self.fname = fname
        self.structure_format = structure_format

        # parse structure data
        self.atoms = None
        if self.structure_format == 'data':
            self.data = DATAReader(fname)
        elif self.structure_format == 'xyz':
            self.data = XYZ2DATAConverter(fname).convert(return_reader=True)

        self.atoms = self.data.atoms
        self.atom_ids = self.atoms.atom_ids
        self.atom_coords = self.atoms.get_coords(as_dict=True)
        self.Nvac = 0
        self.vac_ids = None

    def _random_vacancy_generator(self):
        """Generate random vacancies in structure data."""
        self.vac_ids = \
            np.random.choice(self.atom_ids, size=self.Nvac, replace=False)

    def generate_vacancy_structure(self, show_vmd_selection_cmd=True):
        """Generate vacancy structure.

        Parameters
        ----------
        show_vmd_selection_cmd : bool, optional
            Generate a VMD selection string that can be used to
            select the atoms surrounding the vacancies.

        """
        removed_atoms = self.atoms.filter_atoms(self.vac_ids, invert=False)

        if show_vmd_selection_cmd:
            vmd_vac_selection_cmd = []
            for atom in removed_atoms:
                selection_cmd = "(((x - {:.6f})^2 + ".format(atom.x) + \
                                "(y - {:.6f})^2 + ".format(atom.y) + \
                                "(z - {:.6f})^2) <= 9)".format(atom.z)
                vmd_vac_selection_cmd.append(selection_cmd)

            vmd_vac_selection_cmd = ' or '.join(vmd_vac_selection_cmd)
            print('copy and paste the following VMD command to select\n'
                  'the atoms surrounding the vacancies:\n\n'
                  '{}\n'.format(vmd_vac_selection_cmd))

        remaining_atoms = self.atoms.filter_atoms(self.vac_ids, invert=True)
        vacancy_structure_fname = \
            os.path.splitext(os.path.basename(self.fname))[0] + \
            '+{}_vacancies'.format(self.Nvac)
        DATAWriter.write(fname=vacancy_structure_fname,
                         atoms=remaining_atoms,
                         boxbounds=self.data.boxbounds,
                         comment_line=self.data.comment_line)
        XYZWriter.write(fname=vacancy_structure_fname,
                        atoms=remaining_atoms,
                        comment_line=self.data.comment_line)


class GrapheneVacancyGenerator(VacancyGenerator):
    """Generate vacancies in graphene.

    Parameters
    ----------
    fname : str, optional
        Structure data filename. If you don't provide a structure data file,
        the you **must** provide the structure data parameters to
        generate the structure data file.
    structure_format : {None, str}, optional
        Chemical file format of saved structure data.
        If ``None``, then guess based on ``fname`` file extension.
        Otherwise, must be one of:

            - xyz
            - data

    width : float
        Width of graphene sheet in **nanometers**
    length : float
        Length of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the ``length`` of the sheet sheet.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2.
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    verbose : bool, optional
        Verbose output

    Examples
    --------

    First, we'll generate a single layer sheet of graphene using
    the :py:class:`~sknano.nanogen.GrapheneGenerator` class and then
    "add" some vacancies to it using the
    :py:class:`~sknano.nanogen.GrapheneVacancyGenerator` class.

    >>> from sknano.nanogen import GrapheneGenerator, GrapheneVacancyGenerator
    >>> ZZgraphene_1layer = GrapheneGenerator(width=5, length=10, edge='ZZ')
    >>> ZZgraphene_1layer.save_data(fname='5nmx10nm_ZZ_1layer.data',
    ...                             structure_format='data')

    Now we create an instance of the
    :py:class:`~sknano.nanogen.GrapheneVacancyGenerator` class
    and use the saved structure data output above as input to
    :py:class:`~sknano.nanogen.GrapheneVacancyGenerator`.

    >>> gvacgen = GrapheneVacancyGenerator(fname='5nmx10nm_ZZ_1layer.data',
    ...                                    structure_format='data')

    Now let's add 30 *random* vacancies to the graphene layer.

    .. code-block:: python

       >>> gvacgen.generate_vacancy_structure(Nvac=30, random=True)
       copy and paste the following VMD command to select
       the atoms surrounding the vacancies:

       (((x - -24.157000)^2 + (y - 0.000000)^2 + (z - 37.533974)^2) <= 9) or
       (((x - -24.157000)^2 + (y - 0.000000)^2 + (z - 22.766509)^2) <= 9) or
       (((x - -22.025500)^2 + (y - 0.000000)^2 + (z - 6.768422)^2) <= 9) or
       (((x - -22.025500)^2 + (y - 0.000000)^2 + (z - -7.999044)^2) <= 9) or
       (((x - -24.867500)^2 + (y - 0.000000)^2 + (z - -12.921532)^2) <= 9) or
       (((x - -17.762500)^2 + (y - 0.000000)^2 + (z - 46.148329)^2) <= 9) or
       (((x - -19.894000)^2 + (y - 0.000000)^2 + (z - 30.150241)^2) <= 9) or
       (((x - -18.473000)^2 + (y - 0.000000)^2 + (z - -43.687085)^2) <= 9) or
       (((x - -15.631000)^2 + (y - 0.000000)^2 + (z - 47.378951)^2) <= 9) or
       (((x - -16.341500)^2 + (y - 0.000000)^2 + (z - 16.613398)^2) <= 9) or
       (((x - -14.210000)^2 + (y - 0.000000)^2 + (z - 15.382776)^2) <= 9) or
       (((x - -14.210000)^2 + (y - 0.000000)^2 + (z - -28.919619)^2) <= 9) or
       (((x - -9.947000)^2 + (y - 0.000000)^2 + (z - 32.611486)^2) <= 9) or
       (((x - -9.236500)^2 + (y - 0.000000)^2 + (z - 26.458375)^2) <= 9) or
       (((x - -11.368000)^2 + (y - 0.000000)^2 + (z - 15.382776)^2) <= 9) or
       (((x - -12.078500)^2 + (y - 0.000000)^2 + (z - 4.307177)^2) <= 9) or
       (((x - -7.815500)^2 + (y - 0.000000)^2 + (z - 19.074643)^2) <= 9) or
       (((x - -4.973500)^2 + (y - 0.000000)^2 + (z - -47.378951)^2) <= 9) or
       (((x - -3.552500)^2 + (y - 0.000000)^2 + (z - 19.074643)^2) <= 9) or
       (((x - -1.421000)^2 + (y - 0.000000)^2 + (z - -16.613398)^2) <= 9) or
       (((x - -2.842000)^2 + (y - 0.000000)^2 + (z - -19.074643)^2) <= 9) or
       (((x - -0.710500)^2 + (y - 0.000000)^2 + (z - -39.995218)^2) <= 9) or
       (((x - 12.078500)^2 + (y - 0.000000)^2 + (z - 6.768422)^2) <= 9) or
       (((x - 15.631000)^2 + (y - 0.000000)^2 + (z - 25.227753)^2) <= 9) or
       (((x - 13.499500)^2 + (y - 0.000000)^2 + (z - 16.613398)^2) <= 9) or
       (((x - 15.631000)^2 + (y - 0.000000)^2 + (z - 12.921532)^2) <= 9) or
       (((x - 14.210000)^2 + (y - 0.000000)^2 + (z - -21.535887)^2) <= 9) or
       (((x - 16.341500)^2 + (y - 0.000000)^2 + (z - -35.072730)^2) <= 9) or
       (((x - 18.473000)^2 + (y - 0.000000)^2 + (z - 30.150241)^2) <= 9) or
       (((x - 19.894000)^2 + (y - 0.000000)^2 + (z - 12.921532)^2) <= 9)

    The rendered structure, with vacancies highlighted, looks like:

    .. image:: /images/5nmx10nm_ZZ_1layer+30_vacancies.png

    """
    def __init__(self, fname=None, structure_format=None,
                 width=None, length=None, edge='armchair',
                 element1='C', element2='C', bond=CCbond,
                 nlayers=1, layer_spacing=3.35, stacking_order='AB',
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None, verbose=False):

        if fname is None and width is not None and length is not None:
            try:
                gg = GrapheneGenerator(width=width, length=length,
                                       edge=edge, element1=element1,
                                       element2=element2, bond=bond,
                                       nlayers=nlayers,
                                       layer_spacing=layer_spacing,
                                       stacking_order=stacking_order,
                                       verbose=verbose)
                gg.save_data(structure_format='data')
                fname = gg.fname
            except GrapheneGeneratorError:
                raise VacancyGeneratorError('invalid parameters')

        super(GrapheneVacancyGenerator, self).__init__(
            fname=fname, structure_format=structure_format, verbose=verbose)

    def generate_vacancy_structure(self, Nvac=int, random=False, uniform=False,
                                   bin_axis=None,
                                   show_vmd_selection_cmd=True):
        """Generate vacancy structure.

        Parameters
        ----------
        Nvac : int
            total number of vacancies to "add" to structure data.
        random : bool, optional
            Generate random vacancies in structure data.
        uniform : bool, optional
            Generate uniform vacancies in structure data.
        bin_axis : {'x', 'y', 'z'}, optional
            axis along which to generate uniform vacancies
        show_vmd_selection_cmd : bool, optional

        """
        self.Nvac = Nvac

        if random:
            super(GrapheneVacancyGenerator, self)._random_vacancy_generator()
        else:
            self.vac_ids = np.empty(0, dtype=int)

            # find the coords of each layer
            y_coords = self.atom_coords['y']
            Nlayers = len(set(y_coords))
            layer_coords = np.asarray(sorted(list(set(y_coords))))

            Nvac_per_layer = int(Nvac / Nlayers)
            extra = Nvac % Nlayers
            nbins = Nlayers * Nvac_per_layer

            bin_sets = []
            for n in xrange(Nlayers):
                bin_sets.append(np.arange(n, nbins, Nlayers))
            bin_set_iter = itertools.cycle((bin_sets))

            vac_bin_edges = np.linspace(self.atom_coords[bin_axis].min(),
                                        self.atom_coords[bin_axis].max(),
                                        num=nbins+1)
            vac_coords_along_bin_axis = \
                vac_bin_edges[:-1] + np.diff(vac_bin_edges) / 2

            if self.verbose:
                print('Nvac: {}'.format(Nvac))
                print('Nlayers: {}'.format(Nlayers))
                print('Nvac_per_layer: {}'.format(Nvac_per_layer))
                print('extra vacancies: {}'.format(extra))
                print('nbins: {}'.format(nbins))
                print('bin_sets:\n{}'.format(bin_sets))
                print('vac_bin_edges:'
                      '\n{}'.format(vac_bin_edges))
                print('vac_coords_along_bin_axis:'
                      '\n{}'.format(vac_coords_along_bin_axis))

            for layer_pos in layer_coords:
                bin_set = bin_set_iter.next()
                for vac_pos in vac_coords_along_bin_axis[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(
                            (self.atom_coords['x'] >=
                             (self.atom_coords['x'].min() + 2.5)) &
                            (self.atom_coords['x'] <=
                             (self.atom_coords['x'].max() - 2.5)) &
                            (np.abs(self.atom_coords['y'] - layer_pos)
                                <= 0.5) &
                            (np.abs(self.atom_coords[bin_axis] - vac_pos)
                                <= 1))
                    candidate_vac_atom_ids = \
                        self.atom_ids[candidate_vac_atom_indices]

                    if self.verbose:
                        print('candidate_vac_atom_ids: '
                              '{}\n'.format(candidate_vac_atom_ids))

                    self.vac_ids = \
                        np.r_[self.vac_ids, np.random.choice(
                            candidate_vac_atom_ids)]

        super(GrapheneVacancyGenerator, self).generate_vacancy_structure(
            show_vmd_selection_cmd=show_vmd_selection_cmd)


class NanotubeVacancyGenerator(VacancyGenerator):
    """Generate vacancies in nanotubes.

    Parameters
    ----------
    fname : str, optional
        Structure data filename. If you don't provide a structure data file,
        the you **must** provide the structure data parameters to
        generate the structure data file.
    structure_format : {None, str}, optional
        Chemical file format of saved structure data.
        If ``None``, then guess based on ``fname`` file extension.
        Otherwise, must be one of:

            - xyz
            - data

    n, m : int, optional
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes
    bundle_packing : {None, 'hexagonal', 'cubic'}, optional
        close packing arrangement of bundles
    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    verbose : bool, optional
        Verbose output

    Examples
    --------

    """
    def __init__(self, fname=None, structure_format=None,
                 n=None, m=None, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 vdw_spacing=3.4, bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, rotate_structure=False,
                 rotation_angle=None, rotation_axis=None, verbose=False):

        self.Ntubes = None
        self.Natoms_per_tubes = None

        if fname is None and n is not None and m is not None:
            try:
                ntbg = NanotubeBundleGenerator(n=n, m=m, nx=nx, ny=ny, nz=nz,
                                               element1=element1,
                                               element2=element2,
                                               bond=bond,
                                               vdw_spacing=vdw_spacing,
                                               bundle_packing=bundle_packing,
                                               bundle_geometry=bundle_geometry,
                                               Lx=Lx, Ly=Ly, Lz=Lz,
                                               verbose=verbose)
                ntbg.save_data(structure_format='data')
                fname = ntbg.fname
                self.Ntubes = ntbg.Ntubes
                self.Natoms_per_tube = ntbg.Natoms_per_tube
            except NanotubeGeneratorError:
                raise VacancyGeneratorError('invalid parameters')

        super(NanotubeVacancyGenerator, self).__init__(
            fname=fname, structure_format=structure_format)

    def generate_vacancy_structure(self, Nvac=int, random=False, uniform=False,
                                   bin_axis=None, Ntubes=None,
                                   show_vmd_selection_cmd=True):
        """Generate vacancy structure.

        Parameters
        ----------
        Nvac : int
            total number of vacancies to "add" to structure data.
        random : bool, optional
            Generate random vacancies in structure data.
        uniform : bool, optional
            Generate uniform vacancies in structure data.
        bin_axis : {'x', 'y', 'z'}, optional
            axis along which to generate uniform vacancies
        Ntubes : {None, int}, optional
        show_vmd_selection_cmd : bool, optional

        """
        self.Nvac = Nvac

        if random:
            super(NanotubeVacancyGenerator, self)._random_vacancy_generator()
        else:
            self.vac_ids = np.empty(0, dtype=int)

            if Ntubes is None and self.Ntubes is None:
                raise VacancyGeneratorError('please specify `Ntubes`')
            elif Ntubes is None:
                Ntubes = self.Ntubes

            Natoms = len(self.atom_ids)
            if self.Natoms_per_tube is not None:
                Natoms_per_tube = self.Natoms_per_tube
            else:
                Natoms_per_tube = int(Natoms / Ntubes)

            Nvac_per_tube = int(Nvac / Ntubes)
            extra = Nvac % Ntubes
            nbins = Ntubes * Nvac_per_tube

            bin_sets = []
            for n in xrange(Ntubes):
                bin_sets.append(np.arange(n, Nvac, Ntubes))
            bin_set_iter = itertools.cycle((bin_sets))

            vac_bin_edges = np.linspace(self.atom_coords[bin_axis].min(),
                                        self.atom_coords[bin_axis].max(),
                                        num=nbins+1)
            vac_coords_along_bin_axis = \
                vac_bin_edges[:-1] + np.diff(vac_bin_edges) / 2

            if self.verbose:
                print('Natoms: {}'.format(Natoms))
                print('Natoms_per_tube: {}'.format(Natoms_per_tube))
                print('Nvac: {}'.format(Nvac))
                print('Ntubes: {}'.format(Ntubes))
                print('Nvac_per_tube: {}'.format(Nvac_per_tube))
                print('extra vacancies: {}'.format(extra))
                print('nbins: {}'.format(nbins))
                print('bin_sets:\n{}'.format(bin_sets))
                print('vac_bin_edges:'
                      '\n{}'.format(vac_bin_edges))
                print('vac_coords_along_bin_axis:'
                      '\n{}'.format(vac_coords_along_bin_axis))

            for i in xrange(Ntubes):
                tube_atom_indices = \
                    np.where((self.atom_ids > (Natoms_per_tube * i)) &
                             (self.atom_ids <= (Natoms_per_tube * (i + 1))))

                tube_atom_ids = self.atom_ids[tube_atom_indices]

                tube_coords = self.atoms.get_filtered_coords(tube_atom_ids,
                                                             as_dict=True,
                                                             invert=False)

                if self.verbose:
                    print('tube n: {:d}'.format(i+1))
                    print('tube_atom_ids:\n{}'.format(tube_atom_ids))
                    print('tube_atom_coords:\n{}'.format(tube_coords))

                bin_set = bin_set_iter.next()
                for vac_pos in vac_coords_along_bin_axis[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(np.abs(tube_coords[bin_axis] - vac_pos) <= 1)
                    candidate_vac_atom_ids = \
                        tube_atom_ids[candidate_vac_atom_indices]

                    if self.verbose:
                        print(u'vac_pos along {}-axis: {:.2f} \u00c5'.format(
                            bin_axis, vac_pos))
                        print('N candidate_vac_atom_ids: '
                              '{}\n'.format(len(candidate_vac_atom_ids)))

                    self.vac_ids = \
                        np.r_[self.vac_ids,
                              np.random.choice(candidate_vac_atom_ids)]

        super(NanotubeVacancyGenerator, self).generate_vacancy_structure(
            show_vmd_selection_cmd=show_vmd_selection_cmd)
