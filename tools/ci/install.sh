#!/usr/bin/env bash
set -ex

PYTHON_VERSION="$(python --version)"
echo ${PYTHON_VERSION}

pip install --upgrade pip setuptools
travis_wait pip install -r requirements.txt
python setup.py install
