.. image:: https://travis-ci.org/androomerrill/scikit-nano.svg?branch=master
   :target: https://travis-ci.org/androomerrill/scikit-nano

.. image:: https://coveralls.io/repos/androomerrill/scikit-nano/badge.svg?branch=master&service=github
  :target: https://coveralls.io/github/androomerrill/scikit-nano?branch=master


===========
scikit-nano
===========

*scikit-nano* is a python toolkit for generating and analyzing
nanostructure data.

*scikit-nano* can generate structure data (i.e., atomic coordinates)
for the following classes of nanostructures:

    * Fullerenes
    * Graphene

        * *N*-layer graphene
        * Bilayer graphene with more fine control over relative layer
          orientation, including relative rotation and stacking arrangements.

    * Nanotubes

        * Single-walled nanotubes (SWNTs)
        * SWNT *bundles*
        * Multi-walled nanotubes (MWNTs)
        * MWNT *bundles*


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


Important links
===============

* Documentation: http://scikit-nano.org/doc
* Source code repo: https://github.com/androomerrill/scikit-nano
* Download releases: https://github.com/androomerrill/scikit-nano/releases
* Issue tracker: https://github.com/androomerrill/scikit-nano/issues
* PyPI page: https://pypi.python.org/pypi/scikit-nano

Dependencies
============

* `Python 3.4+ <http://python.org/download/>`_
* `numpy 1.8+ <http://sourceforge.net/projects/numpy/files/NumPy/>`_
* `scipy 0.13+ <http://sourceforge.net/projects/scipy/files/scipy/>`_


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
