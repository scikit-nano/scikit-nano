# -*- coding: utf-8 -*-
"""
===============================================================================
Generate vacancies in structure data (:mod:`sknano.nanogen._vacancy_generator`)
===============================================================================

Module for generating vacancies in nano-structures.

.. currentmodule:: sknano.nanogen._vacancy_generator

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

import itertools
import os

import numpy as np

from ..chemistry import Atoms
from ..nanogen import GrapheneGenerator, GrapheneGeneratorError, \
    NanotubeBundleGenerator, NanotubeGeneratorError
from ..structure_io import DATAReader, DATAWriter, XYZWriter, \
    XYZ2DATAConverter, StructureFormatError, supported_structure_formats
from ..tools.refdata import CCbond

vac_type_cluster_size_map = {'double': 2, 'triple': 3}

__all__ = ['GrapheneVacancyGenerator',
           'NanotubeVacancyGenerator',
           'VacancyGenerator',
           'VacancyGeneratorError']


class VacancyGeneratorError(Exception):
    """Base class for VacancyGenerator exceptions."""
    pass


class VacancyGenerator(object):
    """Base class for generating vacancies in structure data.

    Parameters
    ----------
    fname : str
        structure data filename
    structure_format : {None, str}, optional
        chemical file format of saved structure data.
        If `None`, then guess based on `fname` file extension.
        Otherwise, must be one of:

            - xyz
            - data

    verbose : bool, optional
        Verbose output

    """
    def __init__(self, fname=str, structure_format=None, verbose=False):

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

        self._fname = fname
        self._structure_format = structure_format
        self._verbose = verbose

        # parse structure data
        if self._structure_format == 'data':
            self._structure_data = DATAReader(fname)
        elif self._structure_format == 'xyz':
            self._structure_data = \
                XYZ2DATAConverter(fname).convert(return_reader=True)

        self._atoms = self._structure_data.atoms
        self._atom_ids = self._atoms.atom_ids
        self._atom_coords = self._atoms.get_coords(as_dict=True)

        self._Nvacs = 0
        self._Nvac_clusters = 0
        self._Nvac_sites = 0

        self._vac_ids = np.empty(0, dtype=int)
        self._vac_type = 'single'
        self._cluster_size = 1

        self._vmd_selection_radius = np.sqrt(12)
        self._show_vmd_selection_cmd = True

    @property
    def atoms(self):
        return self._atoms

    @property
    def atom_ids(self):
        return self._atom_ids

    @property
    def Nvacs(self):
        return self._Nvacs

    @property
    def Nvac_clusters(self):
        return self._Nvac_clusters

    @Nvac_clusters.setter
    def Nvac_clusters(self, value):
        self._Nvac_clusters = value

    @property
    def Nvac_sites(self):
        return self._Nvac_sites

    @Nvac_sites.setter
    def Nvac_sites(self, value):
        self._Nvac_sites = value

    @property
    def cluster_size(self):
        return self._cluster_size

    @cluster_size.setter
    def cluster_size(self, value):
        self._cluster_size = value

    @property
    def vac_type(self):
        return self._vac_type

    @vac_type.setter
    def vac_type(self, value):
        self._vac_type = value

    @property
    def vmd_selection_radius(self):
        return self._vmd_selection_radius

    @vmd_selection_radius.setter
    def vmd_selection_radius(self, value):
        self._vmd_selection_radius = value

    @property
    def show_vmd_selection_cmd(self):
        return self._show_vmd_selection_cmd

    @show_vmd_selection_cmd.setter
    def show_vmd_selection_cmd(self, value):
        self._show_vmd_selection_cmd = value

    def _random_vacancy_generator(self):
        """Generate random vacancies in structure data."""
        self._vac_ids = \
            np.random.choice(self._atom_ids,
                             size=self._Nvac_sites,
                             replace=False)

    def _generate_vmd_selection_cmd(self):

        selection_radius = self._vmd_selection_radius
        selections = []
        for atom in self._removed_atoms:
            selection_cmd = \
                "(((x - {:.6f})^2 + ".format(atom.x) + \
                "(y - {:.6f})^2 + ".format(atom.y) + \
                "(z - {:.6f})^2) <= {:.2f})".format(atom.z,
                                                    selection_radius**2)
            selections.append(selection_cmd)

        vmd_selection_cmd = ' or '.join(selections)
        print('copy and paste the following VMD command to select\n'
              'the atoms surrounding the vacancies:\n\n'
              '{}\n'.format(vmd_selection_cmd))

    def _generate_single_vacancies(self):
        self._removed_atoms = \
            self._atoms.filter_atoms(self._vac_ids, invert=False)

    def _generate_multi_vacancies(self):
        vac_type_properties = {'double': {'cluster_size': 2,
                                          'NN_cutoff': 1.5},
                               'triple': {'cluster_size': 3,
                                          'NN_cutoff': 1.5}}
        vac_props = vac_type_properties[self._vac_type]
        self._cluster_size = vac_props['cluster_size']
        self._atoms.NN_cutoff = vac_props['NN_cutoff']
        self._atoms.update_nearest_neighbors()

        vac_atoms = Atoms()
        for vac_id in self._vac_ids:
            vac_atom = self._atoms.get_atom(atomID=vac_id)
            vac_atoms.append(vac_atom)
            vac_atoms.extend(np.random.choice(vac_atom.NN,
                                              size=self._cluster_size-1,
                                              replace=False).tolist())
        self._removed_atoms = \
            self._atoms.filter_atoms(vac_atoms.atom_ids, invert=False)

    def _generate_stone_wales_vacancies(self):
        pass

    def _generate_vacancy_structure(self):
        """Generate vacancy structure."""
        if self._vac_type in ('double', 'triple'):
            self._generate_multi_vacancies()
        elif self._vac_type in ('stone-wales', 'SW'):
            self._generate_stone_wales_vacancies()
        else:
            self._generate_single_vacancies()

        if self._show_vmd_selection_cmd:
            self._generate_vmd_selection_cmd()

        self._Nvacs = self._removed_atoms.Natoms

        self._remaining_atoms = \
            self._atoms.filter_atoms(self._removed_atoms.atom_ids, invert=True)
        #remaining_atoms.assign_unique_ids()

        self._save_vacancy_structure_data()

    def _generate_output_fname(self):
        self._output_fname = \
            os.path.splitext(os.path.basename(self._fname))[0] + \
            '+{}_vacancies'.format(self._Nvacs)

    def _save_vacancy_structure_data(self):
        self._generate_output_fname()
        DATAWriter.write(fname=self._output_fname,
                         atoms=self._remaining_atoms,
                         boxbounds=self._structure_data.boxbounds,
                         comment_line=self._structure_data.comment_line)
        XYZWriter.write(fname=self._output_fname,
                        atoms=self._remaining_atoms,
                        comment_line=self._structure_data.comment_line)


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
        If `None`, then guess based on `fname` file extension.
        Otherwise, must be one of:

            - xyz
            - data

    width : float, optional
        Width of graphene sheet in **nanometers**
    length : float, optional
        Length of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the `length` of the sheet.
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

    Raises
    ------
    :py:class:`VacancyGeneratorError`

    Examples
    --------

    Import the :py:class:`GrapheneVacancyGenerator` class.

    >>> from sknano.nanogen import GrapheneVacancyGenerator

    You can supply any existing structure data file
    (as long as its a supported format) to
    :py:class:`GrapheneVacancyGenerator` using the
    `fname` keyword argument. For example, if you have an existing
    LAMMPS structure data file for graphene named
    *5nmx10nm_ZZ_1layer.data*, you could load that into your
    vacancy generator by typing:

    >>> gvacgen = GrapheneVacancyGenerator(fname='5nmx10nm_ZZ_1layer.data')

    If say, you had your LAMMPS data file saved with a `.lammps`
    extension, then you would need to specify the `structure_format`
    keyword argument as well:

    >>> gvacgen = GrapheneVacancyGenerator(fname='5nmx10nm_ZZ_1layer.lammps',
    ...                                    structure_format='data')

    If you don't have an existing structure data file and want to generate
    your graphene structure on-the-fly, then you need to provide
    structure parameters to the :py:class:`GrapheneVacancyGenerator`
    constructor. These arguments are passed to an instance of
    :py:class:`GrapheneGenerator` "under the hood". See the
    :py:class:`GrapheneGenerator` docs for examples and more
    detailed documentation. The `width` and `length` parameters are
    required, while all others supply default values. In the following
    example, we'll use the :py:class:`GrapheneVacancyGenerator` class
    to generate a **5nm x 10nm**, single layer graphene with `zigzag` edges,
    and then poke some holes in it.

    >>> gvacgen = GrapheneVacancyGenerator(width=5, length=10, edge='ZZ')

    Now let's add 30 *random* vacancies to the graphene layer by
    calling the :py:meth:`~GrapheneVacancyGenerator.generate_vacancy_structure`
    method:

    .. code-block:: python

       >>> gvacgen.generate_vacancy_structure(Nvac_sites=30, random=True)
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

    This will save the vacancy structure data in both
    LAMMPS `data` and `xyz` formats. You can load the `xyz` file
    into VMD and then copy the VMD selection command from the output above
    and paste it into VMD to highlight the atoms surrounding the vacancies.

    The rendered structure, with vacancies highlighted, looks like:

    .. image:: /images/5nmx10nm_ZZ_1layer+30vacancies.png

    """
    def __init__(self, fname=None, structure_format=None,
                 width=None, length=None, edge=None,
                 element1='C', element2='C', bond=CCbond,
                 nlayers=1, layer_spacing=3.35, stacking_order='AB',
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None, verbose=False):

        if fname is None and width is not None and length is not None:
            try:
                gg = GrapheneGenerator(width=width, length=length, edge=edge,
                                       element1=element1, element2=element2,
                                       bond=bond, nlayers=nlayers,
                                       layer_spacing=layer_spacing,
                                       stacking_order=stacking_order,
                                       verbose=verbose)
                gg.save_data(structure_format='data')
                fname = gg.fname
            except GrapheneGeneratorError:
                raise VacancyGeneratorError('invalid parameters')

        super(GrapheneVacancyGenerator, self).__init__(
            fname=fname, structure_format=structure_format, verbose=verbose)

    def generate_vacancy_structure(self, Nvac_sites=None, Nvac_clusters=None,
                                   cluster_size=None, vac_type='single',
                                   uniform=False, bin_axis='z',
                                   distribute_evenly=False,
                                   vmd_selection_radius=np.sqrt(12),
                                   show_vmd_selection_cmd=True):
        """Generate vacancy structure.

        Parameters
        ----------
        Nvac_sites : int, optional
            total number of vacancy sites to "add" to structure data.
        Nvac_clusters : int, optional
            total number of vacancy cluster sites to "add" to structure data.
        uniform : bool, optional
            Generate vacancies uniformly distributed along `bin_axis`.
        show_vmd_selection_cmd : bool, optional
            Generate a VMD selection string that can be used to
            select the atoms surrounding the vacancies.
        vmd_selection_radius : float, optional
            Cutoff radius for VMD selection command in units of **Angstroms**.

        """
        if Nvac_sites is None and Nvac_clusters is None:
            raise ValueError('`Nvac_sites` or `Nvac_clusters` must '
                             'be an integer.')
        elif Nvac_sites is None and Nvac_clusters is not None:
            Nvac_sites = Nvac_clusters
            if cluster_size is None:
                cluster_size = 1

        self._Nvac_sites = Nvac_sites
        self._Nvac_clusters = Nvac_clusters
        self._cluster_size = cluster_size

        self._vac_type = vac_type
        self._vmd_selection_radius = vmd_selection_radius
        self._show_vmd_selection_cmd = show_vmd_selection_cmd

        if uniform:
            # find the coords of each layer
            y_coords = self._atom_coords['y']
            Nlayers = len(set(y_coords))
            layer_coords = np.asarray(sorted(list(set(y_coords))))

            Nvac_sites_per_layer = int(Nvac_sites / Nlayers)
            extra = Nvac_sites % Nlayers

            nbins = Nvac_sites_per_layer * Nlayers

            bin_sets = []
            for n in xrange(Nlayers):
                bin_sets.append(np.arange(n, nbins, Nlayers))
            bin_set_iter = itertools.cycle((bin_sets))

            bin_edges = np.linspace(self._atom_coords[bin_axis].min(),
                                    self._atom_coords[bin_axis].max(),
                                    num=nbins+1)

            bin_mid_pts = bin_edges[:-1] + np.diff(bin_edges) / 2

            if self._verbose:
                print('Nlayers: {}'.format(Nlayers))
                print('Nvac_sites: {}'.format(Nvac_sites))
                print('Nvac_sites_per_layer: {}'.format(Nvac_sites_per_layer))
                print('extra vacancies: {}'.format(extra))
                print('nbins: {}'.format(nbins))
                print('bin_sets:\n{}'.format(bin_sets))
                print('bin_edges:\n{}'.format(bin_edges))
                print('bin_mid_pts:\n{}'.format(bin_mid_pts))

            ortho_axis = 'x' if bin_axis == 'z' else 'z'
            ortho_mid_pt = self._atom_coords[ortho_axis].min() + \
                (self._atom_coords[ortho_axis].max() -
                 self._atom_coords[ortho_axis].min()) / 2
            atoms_along_ortho_axis = {'+': [], '-': []}

            for y in layer_coords:
                bin_set = bin_set_iter.next()
                for vac_pos in bin_mid_pts[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(
                            (self._atom_coords[ortho_axis] >=
                             (self._atom_coords[ortho_axis].min() + 2.5)) &
                            (self._atom_coords[ortho_axis] <=
                             (self._atom_coords[ortho_axis].max() - 2.5)) &
                            (np.abs(self._atom_coords['y'] - y) <= 0.5) &
                            (np.abs(self._atom_coords[bin_axis] - vac_pos)
                                <= 1))
                    candidate_vac_atom_ids = \
                        self._atom_ids[candidate_vac_atom_indices]

                    if self._verbose:
                        print('candidate_vac_atom_ids: '
                              '{}\n'.format(candidate_vac_atom_ids))

                    rand_vac_atom_id = None
                    while True:
                        rand_vac_atom_id = \
                            np.random.choice(candidate_vac_atom_ids)
                        if distribute_evenly:
                            rand_vac_atom = \
                                self._atoms.get_atoms(asarray=True)[
                                    self._atom_ids == rand_vac_atom_id][0]
                            ortho_pos = getattr(rand_vac_atom, ortho_axis)
                            if ortho_pos >= ortho_mid_pt and \
                                    len(atoms_along_ortho_axis['+']) < \
                                    Nvac_sites_per_layer / 2:
                                atoms_along_ortho_axis['+'].append(
                                    rand_vac_atom)
                                break
                            elif ortho_pos < ortho_mid_pt and \
                                    len(atoms_along_ortho_axis['-']) < \
                                    Nvac_sites_per_layer / 2:
                                atoms_along_ortho_axis['-'].append(
                                    rand_vac_atom)
                                break
                            else:
                                continue
                        else:
                            break

                    self._vac_ids = np.r_[self._vac_ids, rand_vac_atom_id]
        else:
            super(GrapheneVacancyGenerator, self)._random_vacancy_generator()

        super(GrapheneVacancyGenerator, self)._generate_vacancy_structure()


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
        If `None`, then guess based on `fname` file extension.
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

    The :py:class:`NanotubeVacancyGenerator` class behaves nearly
    identically to the :py:class:`GrapheneVacancyGenerator` class.
    You can provide any existing structure data file in any supported
    format or else provide structure parameters to generate a structure
    in real-time. The :py:class:`NanotubeVacancyGenerator` uses
    the :py:class:`NanotubeBundleGenerator` class to generate
    nanotubes. See the structure generator classes provided by the
    :py:mod:`~sknano.nanogen` module for more detailed docs.

    In the next example, we'll generate a nanotube bundle and then
    poke some holes in it.

    >>> from sknano.nanogen import NanotubeVacancyGenerator
    >>> ntvg = NanotubeVacancyGenerator(n=10, m=5, Lz=10, fix_Lz=True,
    ...                                 bundle_geometry='hexagon')

    Notice that I used the `Lz` keyword argument to specify the
    **length** of the nanotube bundle in **nanometers**, and I also
    set `fix_Lz=True`. If you don't set `fix_Lz=True`, then the
    length of the nanotube will get truncated such that it's an integer
    multiple of the unit cell length, which,
    depending on the chirality, may end up being very different than
    the desired length set by `Lz`.

    Next, I add 35 vacancies distributed uniformly along the `z`
    axis. When generating vacancies *uniformly*, the code is written such
    it will divide the total number of vacancies
    evenly among all nanotubes and distribute those vacancies uniformly
    along the specified `bin_axis`. The circumferential location of
    vacancies is selected randomly. Future revisions will allow more precise
    control of the positioning.

    .. code-block:: python

       >>> ntvg.generate_vacancy_structure(Nvac_sites=35, uniform=True,
       ...                                 bin_axis='z')
       copy and paste the following VMD command to select
       the atoms surrounding the vacancies:

       (((x - -9.095133)^2 + (y - -2.248367)^2 + (z - -24.706026)^2) <= 9) or
       (((x - -12.610821)^2 + (y - 5.052035)^2 + (z - -13.964275)^2) <= 9) or
       (((x - -17.956207)^2 + (y - 3.045878)^2 + (z - -5.102331)^2) <= 9) or
       (((x - -18.945873)^2 + (y - 0.000000)^2 + (z - 15.575538)^2) <= 9) or
       (((x - -8.581958)^2 + (y - 0.000000)^2 + (z - 6.176506)^2) <= 9) or
       (((x - -2.854725)^2 + (y - 4.324723)^2 + (z - -22.289132)^2) <= 9) or
       (((x - -3.230898)^2 + (y - 4.051418)^2 + (z - -12.353013)^2) <= 9) or
       (((x - -2.854725)^2 + (y - 4.324723)^2 + (z - -3.491069)^2) <= 9) or
       (((x - -1.153094)^2 + (y - -5.052035)^2 + (z - 6.713594)^2) <= 9) or
       (((x - -4.192292)^2 + (y - -3.045878)^2 + (z - 17.455344)^2) <= 9) or
       (((x - 15.800554)^2 + (y - -4.764954)^2 + (z - -21.483501)^2) <= 9) or
       (((x - 11.308348)^2 + (y - -4.563209)^2 + (z - -11.815925)^2) <= 9) or
       (((x - 16.618640)^2 + (y - 4.324723)^2 + (z - -1.074175)^2) <= 9) or
       (((x - 18.615445)^2 + (y - 1.820809)^2 + (z - 8.324857)^2) <= 9) or
       (((x - 9.315470)^2 + (y - -2.657822)^2 + (z - 18.798063)^2) <= 9) or
       (((x - -1.886743)^2 + (y - -13.298492)^2 + (z - -19.872238)^2) <= 9) or
       (((x - -10.112855)^2 + (y - -15.971318)^2 + (z - -9.130488)^2) <= 9) or
       (((x - -8.483271)^2 + (y - -16.848235)^2 + (z - 0.537088)^2) <= 9) or
       (((x - -11.074249)^2 + (y - -14.965779)^2 + (z - 9.936119)^2) <= 9) or
       (((x - -11.330403)^2 + (y - -14.577723)^2 + (z - 18.798063)^2) <= 9) or
       (((x - -10.784342)^2 + (y - 15.329311)^2 + (z - -18.798063)^2) <= 9) or
       (((x - -6.649470)^2 + (y - 6.743161)^2 + (z - -8.593400)^2) <= 9) or
       (((x - -12.063915)^2 + (y - 11.919900)^2 + (z - 0.537088)^2) <= 9) or
       (((x - -2.213175)^2 + (y - 14.168268)^2 + (z - 10.741750)^2) <= 9) or
       (((x - -12.043054)^2 + (y - 12.384407)^2 + (z - 20.677869)^2) <= 9) or
       (((x - 6.186367)^2 + (y - -17.054960)^2 + (z - -17.186801)^2) <= 9) or
       (((x - 7.577549)^2 + (y - -17.054960)^2 + (z - -6.176506)^2) <= 9) or
       (((x - 11.877172)^2 + (y - -10.541309)^2 + (z - 2.148350)^2) <= 9) or
       (((x - 8.918597)^2 + (y - -16.684854)^2 + (z - 12.353013)^2) <= 9) or
       (((x - 11.074249)^2 + (y - -14.965779)^2 + (z - 23.094763)^2) <= 9) or
       (((x - 2.979573)^2 + (y - 15.329311)^2 + (z - -15.038450)^2) <= 9) or
       (((x - 3.300900)^2 + (y - 8.174408)^2 + (z - 3.491069)^2) <= 9) or
       (((x - 11.877172)^2 + (y - 13.298492)^2 + (z - -5.370875)^2) <= 9) or
       (((x - 4.845319)^2 + (y - 7.154946)^2 + (z - 13.158644)^2) <= 9) or
       (((x - 11.877172)^2 + (y - 13.298492)^2 + (z - 24.706026)^2) <= 9)

    Running that command left me with the saved vacancy structure data
    which I loaded into VMD and used the VMD selection command from above
    to highlight the atoms surrounding the vacancies.

    The rendered structure looks like:

    .. image:: /images/1005_hcp_7tube_hexagon+35_vacancies-001.png

    """
    def __init__(self, fname=None, structure_format=None,
                 n=None, m=None, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 vdw_spacing=3.4, bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None, verbose=False):

        self._Ntubes = None
        self._Natoms_per_tube = None

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
                                               fix_Lz=fix_Lz,
                                               verbose=verbose)
                ntbg.save_data(structure_format='data')
                fname = ntbg.fname
                self._Ntubes = ntbg.Ntubes
                self._Natoms_per_tube = ntbg.Natoms_per_tube
            except NanotubeGeneratorError:
                raise VacancyGeneratorError('invalid parameters')

        super(NanotubeVacancyGenerator, self).__init__(
            fname=fname, structure_format=structure_format)

    def generate_vacancy_structure(self, Nvac_sites=None, Nvac_clusters=None,
                                   cluster_size=None, vac_type='single',
                                   uniform=False, bin_axis='z', Ntubes=None,
                                   vmd_selection_radius=np.sqrt(12),
                                   show_vmd_selection_cmd=True):
        """Generate vacancy structure.

        Parameters
        ----------
        Nvac_sites : int, optional
            total number of vacancy sites to "add" to structure data.
        Nvac_clusters : int, optional
            total number of vacancy cluster sites to "add" to structure data.
        random : bool, optional
            Generate random vacancies in structure data.
        uniform : bool, optional
            Generate uniform vacancies in structure data.
        bin_axis : {'x', 'y', 'z'}, optional
            axis along which to generate uniform vacancies
        Ntubes : {None, int}, optional
        show_vmd_selection_cmd : bool, optional
            Generate a VMD selection string that can be used to
            select the atoms surrounding the vacancies.
        vmd_selection_radius : float, optional
            Cutoff radius for VMD selection command in units of **Angstroms**.

        """

        if Nvac_sites is None and Nvac_clusters is None:
            raise ValueError('`Nvac_sites` or `Nvac_clusters` must be an '
                             'integer.')
        elif Nvac_sites is None and Nvac_clusters is not None:
            Nvac_sites = Nvac_clusters
            if cluster_size is None:
                cluster_size = 1

        self._Nvac_sites = Nvac_sites
        self._Nvac_clusters = Nvac_clusters
        self._cluster_size = cluster_size
        self._vac_type = vac_type
        self._vmd_selection_radius = vmd_selection_radius
        self._show_vmd_selection_cmd = show_vmd_selection_cmd

        if uniform:
            if Ntubes is None and self._Ntubes is None:
                raise VacancyGeneratorError('please specify `Ntubes`')
            elif Ntubes is None:
                Ntubes = self._Ntubes

            Natoms = len(self._atom_ids)
            if self._Natoms_per_tube is not None:
                Natoms_per_tube = self._Natoms_per_tube
            else:
                Natoms_per_tube = int(Natoms / Ntubes)

            Nvac_sites_per_tube = int(Nvac_sites / Ntubes)
            extra = Nvac_sites % Ntubes
            nbins = Ntubes * Nvac_sites_per_tube

            bin_sets = []
            for n in xrange(Ntubes):
                bin_sets.append(np.arange(n, Nvac_sites, Ntubes))
            bin_set_iter = itertools.cycle((bin_sets))

            vac_bin_edges = np.linspace(self._atom_coords[bin_axis].min(),
                                        self._atom_coords[bin_axis].max(),
                                        num=nbins+1)
            vac_coords_along_bin_axis = \
                vac_bin_edges[:-1] + np.diff(vac_bin_edges) / 2

            if self._verbose:
                print('Ntubes: {}'.format(Ntubes))
                print('Natoms: {}'.format(Natoms))
                print('Natoms_per_tube: {}'.format(Natoms_per_tube))
                print('Nvac_sites: {}'.format(Nvac_sites))
                print('Nvac_sites_per_tube: {}'.format(Nvac_sites_per_tube))
                print('extra vacancies: {}'.format(extra))
                print('nbins: {}'.format(nbins))
                print('bin_sets:\n{}'.format(bin_sets))
                print('vac_bin_edges:'
                      '\n{}'.format(vac_bin_edges))
                print('vac_coords_along_bin_axis:'
                      '\n{}'.format(vac_coords_along_bin_axis))

            for i in xrange(Ntubes):
                tube_atom_indices = \
                    np.where((self._atom_ids > (Natoms_per_tube * i)) &
                             (self._atom_ids <= (Natoms_per_tube * (i + 1))))

                tube_atom_ids = self._atom_ids[tube_atom_indices]

                tube_coords = self._atoms.get_filtered_coords(tube_atom_ids,
                                                              as_dict=True,
                                                              invert=False)

                if self._verbose:
                    print('tube n: {:d}'.format(i + 1))
                    print('tube_atom_ids:\n{}'.format(tube_atom_ids))
                    print('tube_atom_coords:\n{}'.format(tube_coords))

                bin_set = bin_set_iter.next()
                for vac_pos in vac_coords_along_bin_axis[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(np.abs(tube_coords[bin_axis] - vac_pos) <= 1)
                    candidate_vac_atom_ids = \
                        tube_atom_ids[candidate_vac_atom_indices]

                    if self._verbose:
                        print(u'vac_pos along {}-axis: {:.2f} \u00c5'.format(
                            bin_axis, vac_pos))
                        print('N candidate_vac_atom_ids: '
                              '{}\n'.format(len(candidate_vac_atom_ids)))

                    self._vac_ids = \
                        np.r_[self._vac_ids,
                              np.random.choice(candidate_vac_atom_ids)]
        else:
            super(NanotubeVacancyGenerator, self)._random_vacancy_generator()

        super(NanotubeVacancyGenerator, self)._generate_vacancy_structure()
