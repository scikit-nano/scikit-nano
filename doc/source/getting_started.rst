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

* Python 2.7+
* TubeGen (for generating graphene and nanotube structures)
* VMD (for visualizing structure data)
* Tachyon (for rendering high quality images)

Installing scikit-nano
----------------------

Assuming you already have a working
`python <http://python.org/download/>`_ installation,
and you've downloaded the scikit code, :command:`cd` into the
scikit-nano folder and then run::

    > python setup.py install

This will install the command-line scripts system wide.

To install the scripts in the user's home folder, try::

    > python setup.py install --user
