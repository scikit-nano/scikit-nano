.. _release-guide:

===========================================
How to make a new release of `scikit-nano`
===========================================

A guide for developers for making a new release.

* Edit :file:`setup.py` and change `ISRELEASED` to `True`.

.. _release-testing:

Testing
=======

* Run all of the regression tests by running ``nosetests``
  at the root of the source tree.

.. _release-branching:

Branching
=========


.. _release-packaging:

Packaging
=========


Update PyPI
===========

This step tells PyPI about the release and uploads a source
tarball. This should only be done with final (non-release-candidate)
releases, since doing so will hide any available stable releases.

You may need to set up your `.pypirc` file as described in the
`distutils register command documentation
<http://docs.python.org/3/distutils/packageindex.html>`_.

Then updating the record on PyPI is as simple as::

    python setup.py register

This will hide any previous releases automatically.

Then, to upload the source tarball::

    rm -rf dist
    python setup.py sdist upload

Documentation updates
=====================



Announcing
==========

Announce the release on scikit-nano-announce, scikit-nano-users, and
scikit-nano-dev.  Final (non-release-candidate) versions should also
be announced on python-announce.  Include a summary of highlights from
the CHANGELOG and/or post the whole CHANGELOG since the last release.
