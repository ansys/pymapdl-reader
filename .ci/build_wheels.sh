#!/bin/bash
# builds python wheels on docker container and tests installation

set -e -x

# build based on python version from args
PYTHON_VERSION="$1"
PYBIN="/opt/python/cp${PYTHON_VERSION//.}-cp${PYTHON_VERSION//.}/bin"

# build, don't install
cd io
"${PYBIN}/pip" install build
"${PYBIN}/python" -m build --wheel
auditwheel repair dist/ansys_mapdl_reader*.whl
rm -f dist/*
mv wheelhouse/*manylinux* dist/
