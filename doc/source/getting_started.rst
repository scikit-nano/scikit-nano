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

I've adopted the *scikit-* naming convention for this project in hopes of
contributing it to the collection of
`SciKit libraries <http://scikits.appspot.com>`_ in the near future.

.. _installation:

Installation
============

Required Dependencies
---------------------

* `Python 2.7+ <http://python.org/download/>`_
* VMD (for visualizing structure data)
* Tachyon (for rendering high quality images)

Installing scikit-nano
----------------------

You can install the latest stable release from the
`Python Package Index <http://pypi.python.org/pypi/sknano>`_
using :command:`pip`::

    > pip install sknano

Alternatively you can download a source code tarball from
http://pypi.python.org/pypi/sknano or clone the source code
from the `github repository <http://github.com/androomerrill/scikit-nano>`_
using :command:`git`::

    > git clone https://github.com/androomerrill/scikit-nano.git

:command:`cd` into the source code directory and run::

    > python setup.py install

These commands will probabily fail if you don't have *admin privileges*.
In that case, try installing to the user base directory.
Using :command:`pip`::

    > pip install --user sknano

Or from source::

    > python setup.py install --user
