

MDAtoms
=========================

.. currentmodule:: sknano.core.atoms

.. autoclass:: MDAtoms
   :show-inheritance:

   
     
   

   
   

   .. rubric:: Attributes Summary

   .. autosummary::
   
      ~MDAtoms.CM
      ~MDAtoms.CN_counts
      ~MDAtoms.M
      ~MDAtoms.NNN
      ~MDAtoms.NNrc
      ~MDAtoms.Natoms
      ~MDAtoms.Ntypes
      ~MDAtoms.POAV1
      ~MDAtoms.POAV2
      ~MDAtoms.POAVR
      ~MDAtoms.Rn_counter
      ~MDAtoms.all_angles
      ~MDAtoms.all_bonds
      ~MDAtoms.all_dihedrals
      ~MDAtoms.all_impropers
      ~MDAtoms.angles
      ~MDAtoms.angles_in_degrees
      ~MDAtoms.atom_ids
      ~MDAtoms.atom_tree
      ~MDAtoms.atomtypes
      ~MDAtoms.bonds
      ~MDAtoms.bounding_box
      ~MDAtoms.bounding_region
      ~MDAtoms.bounding_sphere
      ~MDAtoms.bounds
      ~MDAtoms.cell
      ~MDAtoms.cell_matrix
      ~MDAtoms.center_of_mass
      ~MDAtoms.centroid
      ~MDAtoms.charges
      ~MDAtoms.com
      ~MDAtoms.coordinates_bounding_box
      ~MDAtoms.coordination_counts
      ~MDAtoms.coordination_number_counts
      ~MDAtoms.coordination_numbers
      ~MDAtoms.coords
      ~MDAtoms.dihedrals
      ~MDAtoms.distances
      ~MDAtoms.dr
      ~MDAtoms.elements
      ~MDAtoms.etotal
      ~MDAtoms.f
      ~MDAtoms.first_neighbors
      ~MDAtoms.fmtstr
      ~MDAtoms.forces
      ~MDAtoms.fx
      ~MDAtoms.fy
      ~MDAtoms.fz
      ~MDAtoms.i
      ~MDAtoms.ids
      ~MDAtoms.images
      ~MDAtoms.impropers
      ~MDAtoms.indices
      ~MDAtoms.inertia_tensor
      ~MDAtoms.ix
      ~MDAtoms.iy
      ~MDAtoms.iz
      ~MDAtoms.kNN
      ~MDAtoms.ke
      ~MDAtoms.kinetic_energies
      ~MDAtoms.lattice
      ~MDAtoms.lattice_region
      ~MDAtoms.masses
      ~MDAtoms.mol_ids
      ~MDAtoms.molecule_ids
      ~MDAtoms.mols
      ~MDAtoms.moment_of_inertia
      ~MDAtoms.nearest_neighbors
      ~MDAtoms.neighbor_cutoffs
      ~MDAtoms.neighbor_distances
      ~MDAtoms.neighbors
      ~MDAtoms.neighbors_analyzed
      ~MDAtoms.nn_adjacency_list
      ~MDAtoms.nn_adjacency_map
      ~MDAtoms.nn_adjacency_matrix
      ~MDAtoms.nn_seed
      ~MDAtoms.nn_vectors
      ~MDAtoms.pbc
      ~MDAtoms.pe
      ~MDAtoms.positions
      ~MDAtoms.potential_energies
      ~MDAtoms.principal_axes
      ~MDAtoms.principal_moments_of_inertia
      ~MDAtoms.q
      ~MDAtoms.r
      ~MDAtoms.radius_of_gyration
      ~MDAtoms.ring_stats
      ~MDAtoms.rings_per_atom
      ~MDAtoms.rs
      ~MDAtoms.second_neighbors
      ~MDAtoms.serials
      ~MDAtoms.symbols
      ~MDAtoms.third_neighbors
      ~MDAtoms.topology_stats
      ~MDAtoms.total_energies
      ~MDAtoms.typemap
      ~MDAtoms.types
      ~MDAtoms.v
      ~MDAtoms.velocities
      ~MDAtoms.vmd_indices
      ~MDAtoms.volume
      ~MDAtoms.vx
      ~MDAtoms.vy
      ~MDAtoms.vz
      ~MDAtoms.x
      ~MDAtoms.xperiodic
      ~MDAtoms.xs
      ~MDAtoms.y
      ~MDAtoms.yperiodic
      ~MDAtoms.ys
      ~MDAtoms.z
      ~MDAtoms.zperiodic
      ~MDAtoms.zs

   
   

   
   

   .. rubric:: Methods Summary

   .. autosummary::
   
      ~MDAtoms.add_atomtype
      ~MDAtoms.add_atomtypes
      ~MDAtoms.add_ring
      ~MDAtoms.add_type
      ~MDAtoms.add_types
      ~MDAtoms.align_principal_axis
      ~MDAtoms.analyze_POAVs
      ~MDAtoms.analyze_network
      ~MDAtoms.append
      ~MDAtoms.assign_unique_ids
      ~MDAtoms.assign_unique_types
      ~MDAtoms.center_CM
      ~MDAtoms.center_center_of_mass
      ~MDAtoms.center_centroid
      ~MDAtoms.center_com
      ~MDAtoms.clear
      ~MDAtoms.clip_bounds
      ~MDAtoms.compute_POAVs
      ~MDAtoms.copy
      ~MDAtoms.count
      ~MDAtoms.count_neighbors
      ~MDAtoms.count_neighbors_in_self
      ~MDAtoms.extend
      ~MDAtoms.filter
      ~MDAtoms.filter_ids
      ~MDAtoms.filter_vmd_indices
      ~MDAtoms.filtered
      ~MDAtoms.filtered_ids
      ~MDAtoms.filtered_vmd_indices
      ~MDAtoms.get_POAV_attr
      ~MDAtoms.get_angle
      ~MDAtoms.get_angles
      ~MDAtoms.get_atom
      ~MDAtoms.get_atom_from_vmd_index
      ~MDAtoms.get_atoms
      ~MDAtoms.get_atomtypes
      ~MDAtoms.get_bond
      ~MDAtoms.get_bonds
      ~MDAtoms.get_coords
      ~MDAtoms.get_dihedral
      ~MDAtoms.get_dihedrals
      ~MDAtoms.get_improper
      ~MDAtoms.get_impropers
      ~MDAtoms.get_nth_nearest_neighbors
      ~MDAtoms.get_types
      ~MDAtoms.get_vmd_selection_string
      ~MDAtoms.getattr
      ~MDAtoms.index
      ~MDAtoms.insert
      ~MDAtoms.mapatomattr
      ~MDAtoms.pop
      ~MDAtoms.query_atom_tree
      ~MDAtoms.query_ball_point
      ~MDAtoms.query_ball_tree
      ~MDAtoms.query_pairs
      ~MDAtoms.remove
      ~MDAtoms.reset_attrs
      ~MDAtoms.reset_poav_atoms_attrs
      ~MDAtoms.reset_ring_atoms_attrs
      ~MDAtoms.reverse
      ~MDAtoms.rezero
      ~MDAtoms.rezero_coords
      ~MDAtoms.rezero_xyz
      ~MDAtoms.rotate
      ~MDAtoms.select
      ~MDAtoms.set_pbc
      ~MDAtoms.sort
      ~MDAtoms.todict
      ~MDAtoms.translate
      ~MDAtoms.unset_pbc
      ~MDAtoms.update_attrs
      ~MDAtoms.update_neighbor_lists
      ~MDAtoms.update_neighbors
      ~MDAtoms.update_ring_stats
      ~MDAtoms.within_region
      ~MDAtoms.wrap_coords

   
   

   
   

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
   .. autoattribute:: etotal
   .. autoattribute:: f
   .. autoattribute:: first_neighbors
   .. autoattribute:: fmtstr
   .. autoattribute:: forces
   .. autoattribute:: fx
   .. autoattribute:: fy
   .. autoattribute:: fz
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
   .. autoattribute:: ke
   .. autoattribute:: kinetic_energies
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
   .. autoattribute:: pe
   .. autoattribute:: positions
   .. autoattribute:: potential_energies
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
   .. autoattribute:: total_energies
   .. autoattribute:: typemap
   .. autoattribute:: types
   .. autoattribute:: v
   .. autoattribute:: velocities
   .. autoattribute:: vmd_indices
   .. autoattribute:: volume
   .. autoattribute:: vx
   .. autoattribute:: vy
   .. autoattribute:: vz
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

   
   