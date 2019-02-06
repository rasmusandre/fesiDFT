#!/bin/bash
set -e -x

# Install a system package required by our library
yum install -y numpy

ls /opt/python

# Compile wheels
for PYBIN in /opt/python/cp27*/bin; do
    # "${PYBIN}/pip" install --no-binary numpy
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

for PYBIN in /opt/python/cp35*/bin; do
    # "${PYBIN}/pip" install --no-binary numpy
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

for PYBIN in /opt/python/cp36*/bin; do
    # "${PYBIN}/pip" install --no-binary numpy
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
ls wheelhouse/*.whl
for whl in wheelhouse/spglib*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

cp wheelhouse/numpy*whl /io/wheelhouse/

# Install packages and test
for PYBIN in /opt/python/cp27*/bin/; do
    "${PYBIN}/pip" install spglib --no-index -f /io/wheelhouse
done

for PYBIN in /opt/python/cp35*/bin/; do
    "${PYBIN}/pip" install spglib --no-index -f /io/wheelhouse
done

for PYBIN in /opt/python/cp36*/bin/; do
    "${PYBIN}/pip" install spglib --no-index -f /io/wheelhouse
done
