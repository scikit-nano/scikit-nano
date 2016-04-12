:tocdepth: 2

==========================
scikit-nano |version| docs
==========================

:Release: |release|
:Date: |today|

.. |logo_svg| image:: _static/scikit-nano_banner.svg

.. |logo_png| image:: _static/scikit-nano_banner.png

.. raw:: html

   <img src="_images/scikit-nano_banner.svg" onerror="this.src='_images/scikit-nano_banner.png'; this.onerror=null;" width="485"/>

.. only:: latex

    .. image:: _static/scikit-nano_banner.pdf

Welcome to the scikit-nano documentation! scikit-nano is a Python package
intended for nanoscience.

.. _user-docs:

User Documentation
==================

.. only:: html

   :doc:`releases/0.4.0-notes`
   ---------------------------

.. only:: latex

   .. toctree::
      :maxdepth: 1

      releases/0.4.0-notes

scikit-nano quick start guide
------------------------------

.. toctree::
   :maxdepth: 1

   overview
   install
   getting_started
   tutorials/index

Core data structures
----------------------

.. toctree::
   :maxdepth: 1

   core/index
   generators/index

Structure data: Supported formats and I/O
-----------------------------------------

.. toctree::
   :maxdepth: 1

   io/index

scikit-nano data analysis and utilities
-----------------------------------------

.. toctree::
   :maxdepth: 1

   apps/index
   scripts/index
   testing/index
   utils/index

scikit-nano project details
-----------------------------

.. toctree::
   :maxdepth: 1

   releases/index
   known_issues
   license

.. _developer-docs:

Contributing
============

.. toctree::
   :maxdepth: 1

   dev/contribute
   dev/release_guide

Developer Documentation
=======================

.. toctree::
   :maxdepth: 1

   changelog

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _versions:

Other versions
===============

.. ifconfig:: 'dev' in release

   |stable| `documentation <http://docs.scikit-nano.org/>`_

.. ifconfig:: 'dev' not in release

   `scikit-nano development documentation <http://docs.scikit-nano.org/dev>`_
