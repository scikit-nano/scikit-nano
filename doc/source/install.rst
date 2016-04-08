============
Installation
============

Requirements
============

* `Python 3.4+ <http://python.org/download/>`_
* `numpy 1.10+ <http://sourceforge.net/projects/numpy/files/NumPy/>`_
* `scipy 0.16+ <http://sourceforge.net/projects/scipy/files/scipy/>`_

Installing scikit-nano
=======================

You can install the latest stable release from the
`Python Package Index <http://pypi.python.org/pypi/scikit-nano>`_
using **pip**::

    > pip install scikit-nano

Alternatively you can download a source code tarball from
http://pypi.python.org/pypi/scikit-nano or clone the source code
from the `github repo <http://github.com/scikit-nano/scikit-nano>`_
using **git**::

    > git clone https://github.com/scikit-nano/scikit-nano.git

**cd** into the source code directory and run::

    > python setup.py install

These commands will probabily fail if you don't have *admin privileges*.
In that case, try installing to the user base directory.
Using **pip**::

    > pip install --user scikit-nano

Or from source::

    > python setup.py install --user
