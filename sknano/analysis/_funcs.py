# -*- coding: utf-8 -*-
"""
=================================================
Helper functions (:mod:`sknano.analysis._funcs`)
=================================================

.. currentmodule:: sknano.analysis._funcs

"""
from __future__ import division, print_function, absolute_import

__all__ = ['get_nearest_neighbor_counts']


def get_nearest_neighbor_counts(adata=None, datafile=None,
                                filter_elements=None, invert=False):
    pass
#    if adata is not None:
#        atom_data = adata
#    elif datafile is not None and isinstance(datafile, str):
#        atom_data = get_atom_data_asarray(datafile)
#
#    #atom_ids = get_atom_ids_asarray(adata=atom_data)
#
#    filtered_coords = \
#        get_filtered_coords(adata=atom_data, filter_elements=vacIDs)
#    NNtree = spatial.KDTree(filtered_coords)
#    print(NNtree.data)
#
#    # compute things like the distance for which 50% of vacancies
#    # lie within a distance r of each other, etc.
#
#    #return None


def get_nearest_neighbor_comp_array_dict(atom_ids):
    pass
#    nearest_neighbor_comp_dict = OrderedDict()
#    vac_coords = \
#         get_filtered_coords(adata=self.atom_data, filter_elements=vacIDs)
#    #vac_coords = self.get_filtered_coords(vacIDs)
#    for comp in self.position_components:
#        vac_distances = self.get_vac_distances(vac_coords, comp)
#        nearest_neighbor_array = vac_distances[np.where(vac_distances == \
#                                                        vac_distances.min())]
#        nearest_neighbor_comp_dict[comp] = nearest_neighbor_array
#
#    return nearest_neighbor_comp_dict
