#!/usr/bin/env bash
# This script is meant to be called by the "script" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

set -e
mkdir -p $TEST_DIR
# We need the setup.cfg for the nose settings
# cp setup.cfg $TEST_DIR
cd $TEST_DIR

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"

# Skip tests that launch GUIs
export SKNANO_SKIP_GUI_TESTS=1

if [[ "$COVERAGE" == "true" ]]; then
    nosetests -vs --with-coverage sknano
else
    nosetests -vs sknano
fi

# make test-doc test-sphinxext
