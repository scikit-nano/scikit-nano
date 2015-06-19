#!/usr/bin/env python
from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('core', parent_package, top_path)
    config.add_subpackage('atoms')
    config.add_subpackage('crystallography')
    config.add_subpackage('geometric_regions')
    config.add_subpackage('math')
    config.add_subpackage('molecules')
    config.add_subpackage('physics')
    config.add_subpackage('refdata')
    # config.add_subpackage('units')
    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
