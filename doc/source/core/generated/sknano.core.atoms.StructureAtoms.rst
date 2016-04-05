

StructureAtoms
================================

.. currentmodule:: sknano.core.atoms

.. autoclass:: StructureAtoms
   :show-inheritance:

   
     
   

   
   

   .. rubric:: Attributes Summary

   .. autosummary::
   
      ~StructureAtoms.CM
      ~StructureAtoms.CN_counts
      ~StructureAtoms.M
      ~StructureAtoms.NNN
      ~StructureAtoms.NNrc
      ~StructureAtoms.Natoms
      ~StructureAtoms.Ntypes
      ~StructureAtoms.POAV1
      ~StructureAtoms.POAV2
      ~StructureAtoms.POAVR
      ~StructureAtoms.Rn_counter
      ~StructureAtoms.all_angles
      ~StructureAtoms.all_bonds
      ~StructureAtoms.all_dihedrals
      ~StructureAtoms.all_impropers
      ~StructureAtoms.angles
      ~StructureAtoms.angles_in_degrees
      ~StructureAtoms.atom_ids
      ~StructureAtoms.atom_tree
      ~StructureAtoms.atomtypes
      ~StructureAtoms.bonds
      ~StructureAtoms.bounding_box
      ~StructureAtoms.bounding_region
      ~StructureAtoms.bounding_sphere
      ~StructureAtoms.bounds
      ~StructureAtoms.cell
      ~StructureAtoms.cell_matrix
      ~StructureAtoms.center_of_mass
      ~StructureAtoms.centroid
      ~StructureAtoms.charges
      ~StructureAtoms.com
      ~StructureAtoms.coordinates_bounding_box
      ~StructureAtoms.coordination_counts
      ~StructureAtoms.coordination_number_counts
      ~StructureAtoms.coordination_numbers
      ~StructureAtoms.coords
      ~StructureAtoms.dihedrals
      ~StructureAtoms.distances
      ~StructureAtoms.dr
      ~StructureAtoms.elements
      ~StructureAtoms.first_neighbors
      ~StructureAtoms.fmtstr
      ~StructureAtoms.i
      ~StructureAtoms.ids
      ~StructureAtoms.images
      ~StructureAtoms.impropers
      ~StructureAtoms.indices
      ~StructureAtoms.inertia_tensor
      ~StructureAtoms.ix
      ~StructureAtoms.iy
      ~StructureAtoms.iz
      ~StructureAtoms.kNN
      ~StructureAtoms.lattice
      ~StructureAtoms.lattice_region
      ~StructureAtoms.masses
      ~StructureAtoms.mol_ids
      ~StructureAtoms.molecule_ids
      ~StructureAtoms.mols
      ~StructureAtoms.moment_of_inertia
      ~StructureAtoms.nearest_neighbors
      ~StructureAtoms.neighbor_cutoffs
      ~StructureAtoms.neighbor_distances
      ~StructureAtoms.neighbors
      ~StructureAtoms.neighbors_analyzed
      ~StructureAtoms.nn_adjacency_list
      ~StructureAtoms.nn_adjacency_map
      ~StructureAtoms.nn_adjacency_matrix
      ~StructureAtoms.nn_seed
      ~StructureAtoms.nn_vectors
      ~StructureAtoms.pbc
      ~StructureAtoms.positions
      ~StructureAtoms.principal_axes
      ~StructureAtoms.principal_moments_of_inertia
      ~StructureAtoms.q
      ~StructureAtoms.r
      ~StructureAtoms.radius_of_gyration
      ~StructureAtoms.ring_stats
      ~StructureAtoms.rings_per_atom
      ~StructureAtoms.rs
      ~StructureAtoms.second_neighbors
      ~StructureAtoms.serials
      ~StructureAtoms.symbols
      ~StructureAtoms.third_neighbors
      ~StructureAtoms.topology_stats
      ~StructureAtoms.typemap
      ~StructureAtoms.types
      ~StructureAtoms.vmd_indices
      ~StructureAtoms.volume
      ~StructureAtoms.x
      ~StructureAtoms.xperiodic
      ~StructureAtoms.xs
      ~StructureAtoms.y
      ~StructureAtoms.yperiodic
      ~StructureAtoms.ys
      ~StructureAtoms.z
      ~StructureAtoms.zperiodic
      ~StructureAtoms.zs

   
   

   
   

   .. rubric:: Methods Summary

   .. autosummary::
   
      ~StructureAtoms.add_atomtype
      ~StructureAtoms.add_atomtypes
      ~StructureAtoms.add_ring
      ~StructureAtoms.add_type
      ~StructureAtoms.add_types
      ~StructureAtoms.align_principal_axis
      ~StructureAtoms.analyze_POAVs
      ~StructureAtoms.analyze_network
      ~StructureAtoms.append
      ~StructureAtoms.assign_unique_ids
      ~StructureAtoms.assign_unique_types
      ~StructureAtoms.center_CM
      ~StructureAtoms.center_center_of_mass
      ~StructureAtoms.center_centroid
      ~StructureAtoms.center_com
      ~StructureAtoms.clear
      ~StructureAtoms.clip_bounds
      ~StructureAtoms.compute_POAVs
      ~StructureAtoms.copy
      ~StructureAtoms.count
      ~StructureAtoms.count_neighbors
      ~StructureAtoms.count_neighbors_in_self
      ~StructureAtoms.extend
      ~StructureAtoms.filter
      ~StructureAtoms.filter_ids
      ~StructureAtoms.filter_vmd_indices
      ~StructureAtoms.filtered
      ~StructureAtoms.filtered_ids
      ~StructureAtoms.filtered_vmd_indices
      ~StructureAtoms.get_POAV_attr
      ~StructureAtoms.get_angle
      ~StructureAtoms.get_angles
      ~StructureAtoms.get_atom
      ~StructureAtoms.get_atom_from_vmd_index
      ~StructureAtoms.get_atoms
      ~StructureAtoms.get_atomtypes
      ~StructureAtoms.get_bond
      ~StructureAtoms.get_bonds
      ~StructureAtoms.get_coords
      ~StructureAtoms.get_dihedral
      ~StructureAtoms.get_dihedrals
      ~StructureAtoms.get_improper
      ~StructureAtoms.get_impropers
      ~StructureAtoms.get_nth_nearest_neighbors
      ~StructureAtoms.get_types
      ~StructureAtoms.get_vmd_selection_string
      ~StructureAtoms.getattr
      ~StructureAtoms.index
      ~StructureAtoms.insert
      ~StructureAtoms.mapatomattr
      ~StructureAtoms.pop
      ~StructureAtoms.query_atom_tree
      ~StructureAtoms.query_ball_point
      ~StructureAtoms.query_ball_tree
      ~StructureAtoms.query_pairs
      ~StructureAtoms.remove
      ~StructureAtoms.reset_attrs
      ~StructureAtoms.reset_poav_atoms_attrs
      ~StructureAtoms.reset_ring_atoms_attrs
      ~StructureAtoms.reverse
      ~StructureAtoms.rezero
      ~StructureAtoms.rezero_coords
      ~StructureAtoms.rezero_xyz
      ~StructureAtoms.rotate
      ~StructureAtoms.select
      ~StructureAtoms.set_pbc
      ~StructureAtoms.sort
      ~StructureAtoms.todict
      ~StructureAtoms.translate
      ~StructureAtoms.unset_pbc
      ~StructureAtoms.update_attrs
      ~StructureAtoms.update_neighbor_lists
      ~StructureAtoms.update_neighbors
      ~StructureAtoms.update_ring_stats
      ~StructureAtoms.within_region
      ~StructureAtoms.wrap_coords

   
   

   
   

   .. rubric:: Attributes Documentation

   
   .. autoattribute:: CM
   .. autoattribute:: CN_counts
   .. autoattribute:: M
   .. autoattribute:: NNN
   .. autoattribute:: NNrc
   .. autoattribute:: Natoms
   .. autoattribute:: Ntypes
   .. autoattribute:: POAV1
   .. autoattribute:: POAV2
   .. autoattribute:: POAVR
   .. autoattribute:: Rn_counter
   .. autoattribute:: all_angles
   .. autoattribute:: all_bonds
   .. autoattribute:: all_dihedrals
   .. autoattribute:: all_impropers
   .. autoattribute:: angles
   .. autoattribute:: angles_in_degrees
   .. autoattribute:: atom_ids
   .. autoattribute:: atom_tree
   .. autoattribute:: atomtypes
   .. autoattribute:: bonds
   .. autoattribute:: bounding_box
   .. autoattribute:: bounding_region
   .. autoattribute:: bounding_sphere
   .. autoattribute:: bounds
   .. autoattribute:: cell
   .. autoattribute:: cell_matrix
   .. autoattribute:: center_of_mass
   .. autoattribute:: centroid
   .. autoattribute:: charges
   .. autoattribute:: com
   .. autoattribute:: coordinates_bounding_box
   .. autoattribute:: coordination_counts
   .. autoattribute:: coordination_number_counts
   .. autoattribute:: coordination_numbers
   .. autoattribute:: coords
   .. autoattribute:: dihedrals
   .. autoattribute:: distances
   .. autoattribute:: dr
   .. autoattribute:: elements
   .. autoattribute:: first_neighbors
   .. autoattribute:: fmtstr
   .. autoattribute:: i
   .. autoattribute:: ids
   .. autoattribute:: images
   .. autoattribute:: impropers
   .. autoattribute:: indices
   .. autoattribute:: inertia_tensor
   .. autoattribute:: ix
   .. autoattribute:: iy
   .. autoattribute:: iz
   .. autoattribute:: kNN
   .. autoattribute:: lattice
   .. autoattribute:: lattice_region
   .. autoattribute:: masses
   .. autoattribute:: mol_ids
   .. autoattribute:: molecule_ids
   .. autoattribute:: mols
   .. autoattribute:: moment_of_inertia
   .. autoattribute:: nearest_neighbors
   .. autoattribute:: neighbor_cutoffs
   .. autoattribute:: neighbor_distances
   .. autoattribute:: neighbors
   .. autoattribute:: neighbors_analyzed
   .. autoattribute:: nn_adjacency_list
   .. autoattribute:: nn_adjacency_map
   .. autoattribute:: nn_adjacency_matrix
   .. autoattribute:: nn_seed
   .. autoattribute:: nn_vectors
   .. autoattribute:: pbc
   .. autoattribute:: positions
   .. autoattribute:: principal_axes
   .. autoattribute:: principal_moments_of_inertia
   .. autoattribute:: q
   .. autoattribute:: r
   .. autoattribute:: radius_of_gyration
   .. autoattribute:: ring_stats
   .. autoattribute:: rings_per_atom
   .. autoattribute:: rs
   .. autoattribute:: second_neighbors
   .. autoattribute:: serials
   .. autoattribute:: symbols
   .. autoattribute:: third_neighbors
   .. autoattribute:: topology_stats
   .. autoattribute:: typemap
   .. autoattribute:: types
   .. autoattribute:: vmd_indices
   .. autoattribute:: volume
   .. autoattribute:: x
   .. autoattribute:: xperiodic
   .. autoattribute:: xs
   .. autoattribute:: y
   .. autoattribute:: yperiodic
   .. autoattribute:: ys
   .. autoattribute:: z
   .. autoattribute:: zperiodic
   .. autoattribute:: zs

   
   

   
   

   .. rubric:: Methods Documentation

   
   .. automethod:: add_atomtype
   .. automethod:: add_atomtypes
   .. automethod:: add_ring
   .. automethod:: add_type
   .. automethod:: add_types
   .. automethod:: align_principal_axis
   .. automethod:: analyze_POAVs
   .. automethod:: analyze_network
   .. automethod:: append
   .. automethod:: assign_unique_ids
   .. automethod:: assign_unique_types
   .. automethod:: center_CM
   .. automethod:: center_center_of_mass
   .. automethod:: center_centroid
   .. automethod:: center_com
   .. automethod:: clear
   .. automethod:: clip_bounds
   .. automethod:: compute_POAVs
   .. automethod:: copy
   .. automethod:: count
   .. automethod:: count_neighbors
   .. automethod:: count_neighbors_in_self
   .. automethod:: extend
   .. automethod:: filter
   .. automethod:: filter_ids
   .. automethod:: filter_vmd_indices
   .. automethod:: filtered
   .. automethod:: filtered_ids
   .. automethod:: filtered_vmd_indices
   .. automethod:: get_POAV_attr
   .. automethod:: get_angle
   .. automethod:: get_angles
   .. automethod:: get_atom
   .. automethod:: get_atom_from_vmd_index
   .. automethod:: get_atoms
   .. automethod:: get_atomtypes
   .. automethod:: get_bond
   .. automethod:: get_bonds
   .. automethod:: get_coords
   .. automethod:: get_dihedral
   .. automethod:: get_dihedrals
   .. automethod:: get_improper
   .. automethod:: get_impropers
   .. automethod:: get_nth_nearest_neighbors
   .. automethod:: get_types
   .. automethod:: get_vmd_selection_string
   .. automethod:: getattr
   .. automethod:: index
   .. automethod:: insert
   .. automethod:: mapatomattr
   .. automethod:: pop
   .. automethod:: query_atom_tree
   .. automethod:: query_ball_point
   .. automethod:: query_ball_tree
   .. automethod:: query_pairs
   .. automethod:: remove
   .. automethod:: reset_attrs
   .. automethod:: reset_poav_atoms_attrs
   .. automethod:: reset_ring_atoms_attrs
   .. automethod:: reverse
   .. automethod:: rezero
   .. automethod:: rezero_coords
   .. automethod:: rezero_xyz
   .. automethod:: rotate
   .. automethod:: select
   .. automethod:: set_pbc
   .. automethod:: sort
   .. automethod:: todict
   .. automethod:: translate
   .. automethod:: unset_pbc
   .. automethod:: update_attrs
   .. automethod:: update_neighbor_lists
   .. automethod:: update_neighbors
   .. automethod:: update_ring_stats
   .. automethod:: within_region
   .. automethod:: wrap_coords

   
   