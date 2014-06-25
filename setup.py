# -*- coding: utf-8 -*-
"""Python toolkit for generating and analyzing nanostructure data"""

from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext'

import os
import sys
import subprocess

from setuptools import find_packages

from numpy.distutils.core import setup

if sys.version_info[:2] < (2, 7):
    raise RuntimeError("Python version 2.7 required.")

DISTNAME = 'scikit-nano'
DESCRIPTION = __doc__
LONG_DESCRIPTION = open('README.rst').read()
AUTHOR = 'Andrew Merrill'
AUTHOR_EMAIL = 'androomerrill@gmail.com'
MAINTAINER = AUTHOR
MAINTAINER_EMAIL = AUTHOR_EMAIL
URL = 'http://projects.geekcode.io/scikit-nano/doc'
DOWNLOAD_URL = 'https://github.com/androomerrill/scikit-nano'
KEYWORDS = 'nano nano-structures nanostructures nanotubes graphene LAMMPS XYZ'
LICENSE = 'BSD 2-Clause'
CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Programming Language :: Python
Programming Language :: Python :: 2.7
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Chemistry
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Visualization
Topic :: Software Development
Topic :: Software Development :: Libraries :: Python Modules

"""

MAJOR = 0
MINOR = 2
MICRO = 25
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


def git_version():
    """Return the GIT version as a string."""
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = \
            subprocess.Popen(cmd,
                             stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def get_version_info():
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of sknano.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('sknano/version.py'):
        # must be a source distribution, use existing version file
        # load it as a separate module to not load sknano/__init__.py
        import imp
        version = imp.load_source('sknano.version', 'sknano/version.py')
        GIT_REVISION = version.git_revision
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='sknano/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SCIKIT-NANO SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


def configuration(parent_package='', top_path=None):
    if os.path.exists('MANIFEST'):
        os.remove('MANIFEST')

    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)

    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('sknano')

    return config


def setup_package():

    # Rewrite the version file everytime
    write_version_py()

    FULLVERSION, GIT_REVISION = get_version_info()

    setup_options = dict(
        name=DISTNAME,
        version=FULLVERSION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url=URL,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        keywords=KEYWORDS,
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Linux", "OS-X", "Unix", "Windows"],
        test_suite='nose.collector',
        configuration=configuration,
        packages=find_packages(exclude=['doc']),
        include_package_data=True,
        exclude_package_data={'': ['*.gif', '*.html', '*.ui']},
        zip_safe=False,
        install_requires=['numpy>=1.8', 'scipy>=0.14', 'pint>=0.5'],
        entry_points={
            'console_scripts': [
                'nanogen = sknano.scripts.nanogen:main',
                'nanogenui = sknano.scripts.nanogenui:main',
            ]
        }
    )

    setup(**setup_options)

if __name__ == '__main__':
    setup_package()
