#!/usr/bin/env bash
set -ex

PYTHON_VERSION="$(python --version)"
echo ${PYTHON_VERSION}

pip install --upgrade pip setuptools
travis_wait pip install -r requirements.txt
pip install nose

export COVERAGE="true"

if [[ "$COVERAGE" == "true" ]]; then
    pip install coverage coveralls
fi

python setup.py install
