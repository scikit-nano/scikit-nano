#!/usr/bin/env python
"""
A setup.py script to use setuptools, which gives egg goodness, etc.
"""
from __future__ import absolute_import

from setuptools import setup
exec(compile(open('setup.py').read(), 'setup.py', 'exec'))
