#!/bin/bash
# builds python wheels on docker container and tests installation

set -e -x

# build based on python version from args
PYTHON_VERSION="$1"
case $PYTHON_VERSION in
3.7)
  PYBIN="/opt/python/cp37-cp37m/bin"
  ;;
3.8)
  PYBIN="/opt/python/cp38-cp38/bin"
  ;;
3.9)
  PYBIN="/opt/python/cp39-cp39/bin"
  ;;
3.10)
  PYBIN="/opt/python/cp310-cp310/bin"
  ;;
3.11)
  PYBIN="/opt/python/cp311-cp311/bin"
  ;;
esac

# build, don't install
cd io
"${PYBIN}/pip" install build
"${PYBIN}/python" -m build --wheel
auditwheel repair dist/ansys_mapdl_reader*.whl
rm -f dist/*
mv wheelhouse/*manylinux* dist/
