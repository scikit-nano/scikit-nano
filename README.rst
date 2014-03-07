
scikit-nano
===========

*scikit-nano* is a python toolkit for generating and analyzing
nanostructure data.

Currently, its primary utility is generating nanostructure data
(i.e. atomic coordinates) for the following nanostructure materials:

    * Graphene:

        * Single layer graphene
        * Bi-layer graphene with layers rotated relative to each other
          by any angle and different layer stacking arrangements
        * *N*-layer graphene

    * Nanotubes:

        * Single-walled nanotubes (SWNTs)
        * SWNT *bundles*
        * Multi-walled nanotubes (MWNTs)
        * MWNT *bundles*

It currently supports saving structure data in the following formats:

    * `LAMMPS data`
    * `xyz`

Secondary to its structure generating functions are its
*structure analysis tools* including:

    * defect/vacancy structure analysis
    * nearest-neighbor analysis
    * ...


Important links
===============

* Documentation: http://projects.geekcode.io/scikit-nano/doc
* Source code repo: https://github.com/androomerrill/scikit-nano
* Download releases: https://github.com/androomerrill/scikit-nano/releases
* Issue tracker: https://github.com/androomerrill/scikit-nano/issues
* PyPI page: https://pypi.python.org/pypi/scikit-nano

Dependencies
============

Required Dependencies
---------------------
* `Python 2.7+ <http://python.org/download/>`_
* `numpy 1.8+ <http://sourceforge.net/projects/numpy/files/NumPy/>`_

Optional Dependencies
---------------------
* `scipy 0.13+ <http://sourceforge.net/projects/scipy/files/scipy/>`_ (for
  nearest-neighbor analysis)
* `pint 0.4+ <https://pypi.python.org/pypi/Pint/>`_ (for physical units)
* `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_ (for visualizing structure data)
* `Tachyon Ray Tracer <http://jedi.ks.uiuc.edu/~johns/raytracer/>`_ (for
  rendering high quality images)

Installation
=============

You can install the latest stable release from the
`Python Package Index <http://pypi.python.org/pypi/scikit-nano>`_
using **pip**::

    > pip install scikit-nano

Alternatively you can download a source code tarball from
http://pypi.python.org/pypi/scikit-nano or clone the source code
from the `github repo <http://github.com/androomerrill/scikit-nano>`_
using **git**::

    > git clone https://github.com/androomerrill/scikit-nano.git

**cd** into the source code directory and run::

    > python setup.py install

These commands will probabily fail if you don't have *admin privileges*.
In that case, try installing to the user base directory.
Using **pip**::

    > pip install --user scikit-nano

Or from source::

    > python setup.py install --user
