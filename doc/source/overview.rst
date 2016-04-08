========
Overview
========

**scikit-nano** is a Python toolkit for nanoscience.

``scikit-nano`` Package
=======================

The main feature of **scikit-nano** compared to other existing software
and toolkits, is its structure data generators.

*scikit-nano* can generate structure data (i.e., atomic coordinates)
for the following classes of nanostructures:

    * Nano Structures

        * Fullerenes
        * Graphene
        * Nanotubes

            * Single-walled nanotubes (SWNTs) and SWNT bundles.
            * Multi-walled nanotubes (MWNTs) and MWNT bundles

    * Crystal Structures

    * Layered Structures

The following structure data formats are supported:

    * `xyz`
    * `LAMMPS data` (*limited support for full format spec.*)
    * `LAMMPS dump` (*limited support for full format spec.*)


Extending input/output capabilities with more structure data formats
such as *pdb*, *json*, *zmatrix*, etc. is queued for development

Secondary to its structure generating functions are its
*structure analysis tools* including:

    * defect/vacancy structure analysis
    * nearest-neighbor analysis
    * POAV analysis
