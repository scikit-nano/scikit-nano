#!/usr/bin/env bash
set -ex

PYTHON_VERSION="$(python --version)"
echo ${PYTHON_VERSION}

export CC=gcc
export CXX=g++

echo 'List files from cached directories'
echo 'pip:'
ls $HOME/.cache/pip
echo 'download'
ls $HOME/download


if [[ "$DISTRIB" == "conda" ]]; then
    # Deactivate the travis-provided virtual environment and setup a
    # conda-based environment instead
    deactivate

    # Use the miniconda installer for faster download / install of conda
    # itself
    pushd .
    cd
    mkdir -p download
    cd download
    echo "Cached in $HOME/download :"
    ls -l
    echo
    if [[ ! -f miniconda.sh ]]
        then
        wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
            -O miniconda.sh
        fi
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    conda info -a
    popd

    # Configure the conda environment and put it in the path using the
    # provided versions
    if [[ "$INSTALL_MKL" == "true" ]]; then
        conda create -n testenv python=$PYTHON_VERSION pip nose \
            numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION numpy scipy \
            cython=$CYTHON_VERSION matplotlib pandas xlsxwriter \
            hdf5 pytables libgfortran mkl
    else
        conda create -n testenv python=$PYTHON_VERSION pip nose \
            numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION cython=$CYTHON_VERSION \
            matplotlib pandas xlsxwriter hdf5 pytables libgfortran
    fi
    source activate testenv

    # Install nose-timer via pip
    pip install nose-timer
    travis_wait pip -v install -r requirements

    # Resolve MKL usage


elif [[ "$DISTRIB" == "ubuntu" ]]; then
    # At the time of writing numpy 1.9.1 is included in the travis
    # virtualenv but we want to used numpy installed through apt-get
    # install.
    deactivate
    # Create a new virtualenv using system site packages for numpy and scipy
    virtualenv --system-site-packages testvenv
    source testvenv/bin/activate
    pip install nose nose-timer
    pip install cython
    travis_wait pip -v install -r requirements
fi

if [[ "$COVERAGE" == "true" ]]; then
    pip install coverage coveralls
fi

if [ ! -d "$CACHED_BUILD_DIR" ]; then
    mkdir -p $CACHED_BUILD_DIR
fi

rsync -av --exclude '.git/' --exclude='testvenv/' \
      $TRAVIS_BUILD_DIR $CACHED_BUILD_DIR

cd $CACHED_BUILD_DIR/scikit-nano

# Build scikit-learn in the install.sh script to collapse the verbose
# build output in the travis output when it succeeds.
python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python setup.py develop

# pip install --upgrade pip setuptools
# travis_wait pip install -r requirements.txt
# pip install nose

# export COVERAGE="true"

# if [[ "$COVERAGE" == "true" ]]; then
#     pip install coverage coveralls
# fi

# python setup.py install
