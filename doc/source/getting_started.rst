.. _getting_started:

===============
Getting started
===============

.. _introduction:

Introduction
============

.. include:: ../../README.rst
   :start-line: 4
   :end-line: -2

This package is now featured in the collection of
`SciKit libraries <http://scikits.appspot.com/scikits>`_

.. _installation:

Installation
============

Required Dependencies
---------------------
* `Python 2.7+ <http://python.org/download/>`_
* `numpy 1.8+ http://sourceforge.net/projects/numpy/files/NumPy/>`_

Optional Dependencies
---------------------
* `scipy 0.13+ http://sourceforge.net/projects/scipy/files/scipy/>`_ (for
  nearest-neighbor analysis)
* `pint 0.4+ <https://pypi.python.org/pypi/Pint/>`_ (for physical units)
* `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_ (for visualizing structure data)
* `Tachyon Ray Tracer <http://jedi.ks.uiuc.edu/~johns/raytracer/>`_ (for
  rendering high quality images)

Downloading/installing scikit-nano
-----------------------------------
The package source is available from a couple locations:

* `Python Package Index <https://pypi.python.org/pypi/scikit-nano>`_ (latest
  stable release)
* `github repo <https://github.com/androomerrill/scikit-nano>`_ (stable and
  development branches)

You can install the latest stable release from the
`Python Package Index <http://pypi.python.org/pypi/scikit-nano>`_
using :command:`pip`::

    > pip install scikit-nano

Alternatively you can download a source code tarball from
http://pypi.python.org/pypi/scikit-nano or clone the source code
from the `github repo <http://github.com/androomerrill/scikit-nano>`_
using :command:`git`::

    > git clone https://github.com/androomerrill/scikit-nano.git

:command:`cd` into the source code directory and run::

    > python setup.py install

These commands will probabily fail if you don't have *admin privileges*.
In that case, try installing to the user base directory.
Using :command:`pip`::

    > pip install --user scikit-nano

Or from source::

    > python setup.py install --user
