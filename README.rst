===========
scikit-nano
===========

*scikit-nano* is a python toolkit for generating nanostructure data and
nanostructure analysis.

Its primary utility is generating nanostructure data (i.e. atomic coordinates)
for the following nanostructure materials:

    * Graphene:

        * Single layer graphene
        * Bi-layer graphene with layers rotated relative to each other
          by any angle and different layer stacking arrangements
        * :math:`N`-layer graphene

    * Nanotubes:

        * Single-walled nanotubes (SWNTs)
        * SWNT *bundles*
        * Multi-walled nanotubes (MWNTs)
        * MWNT *bundles*

It currently supports saving structure data in the following formats:

    * `xyz`
    * LAMMPS `data`

Secondary to its structure generating functions are its
*structure analysis tools* including:

    * defect/vacancy structure analysis
    * nearest-neighbor analysis
    * ...

I welcome feedback and contributions to development!

For documentation, see:
`scikit-nano documentation <http://projects.geekspin.net/scikit-nano/doc>`_
