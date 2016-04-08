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

    .. image:: _static/scikit-nano_logo.pdf

Welcome to the scikit-nano documentation! scikit-nano is a Python package
intended for nanoscience.

.. _user-docs:

User Documentation
==================

.. only:: html

   :doc:`releases/0.4.0-notes`
   ---------------------

.. only:: latex

   .. toctree::
      :maxdepth: 1

      releases/0.4.0-notes

.. toctree::
   :maxdepth: 1

   install
   getting_started

.. toctree::
   :maxdepth: 2

   tutorial/index

**scikit-nano project details**

.. toctree::
   :maxdepth: 1

   releases/index
   known_issues
   license

Contributing
============

.. toctree::
   :maxdepth: 1

   contribute

.. _developer-docs:

Developer Documentation
=======================

.. toctree::
   :maxdepth: 1

   changelog

.. _reference:

API Reference
=============

.. toctree::
   :maxdepth: 1

   apps/index
   core/index
   generators/index
   io/index
   scripts/index
   testing/index
   utils/index

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
