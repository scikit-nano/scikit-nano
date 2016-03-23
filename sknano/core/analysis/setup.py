#!/usr/bin/env python
from __future__ import division, print_function, absolute_import, \
    unicode_literals
import os


def configuration(parent_package='core', top_path=None):
    import numpy as np
    from numpy.distutils.misc_util import Configuration
    from distutils.sysconfig import get_python_inc
    config = Configuration('analysis', parent_package, top_path)
    config.add_data_dir('tests')

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    include_dirs = [get_python_inc()]
    if include_dirs[0] != get_python_inc(plat_specific=1):
        include_dirs.append(get_python_inc(plat_specific=1))
    include_dirs.append(np.get_include())

    extra_compile_args = []
    #extra_compile_args.append('-std=c++11')

    config.add_extension('_ring_finder',
                         sources=['_ring_finder.cxx'],
                         include_dirs=include_dirs,
                         extra_compile_args=extra_compile_args,
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
