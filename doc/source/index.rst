=========================
scikit-nano documentation
=========================

:Release: |release|
:Date: |today|

.. include:: ../../README.rst
   :start-line: 10
   :end-line: 12

.. toctree::
   :maxdepth: 1

   getting_started.rst

.. toctree::
   :maxdepth: 2

   tutorial/index.rst

.. toctree::
   :maxdepth: 1

   contribute.rst
   credits.rst
   release.rst


.. _reference:

=========
Reference
=========

.. toctree::
   :maxdepth: 1

   apps.rst
   core.rst
   generators.rst
   io.rst
   scripts.rst
   structures.rst
   testing.rst
   utils.rst


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _versions:

Other versions
===============

.. ifconfig:: 'dev' in release

   |stable| `documentation <http://scikit-nano.org/doc>`_

.. ifconfig:: 'dev' not in release

   `scikit-nano development documentation <http://scikit-nano.org/doc/dev>`_
