# Copyright (C) 2021 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import pathlib

from conftest import skip_no_graphics
import numpy as np
import pytest

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples
from ansys.mapdl.reader.full import FullFile
from ansys.mapdl.reader.misc.checks import run_if_scipy_required

try:
    run_if_scipy_required()
    import scipy

    is_scipy_available = True
except ImportError:
    is_scipy_available = False

test_path = os.path.dirname(os.path.abspath(__file__))
testfiles_path = os.path.join(test_path, "testfiles")

skip_scipy = pytest.mark.skipif(not is_scipy_available, reason="Requires SciPy")


@pytest.fixture()
def sparse_full_pathlib_full_file():
    filename = os.path.join(testfiles_path, "sparse.full")
    return FullFile(pathlib.Path(filename))


@skip_no_graphics
@pytest.fixture()
def sparse_full():
    filename = os.path.join(testfiles_path, "sparse.full")
    return pymapdl_reader.read_binary(filename)


@skip_no_graphics
def test_fullreader():
    fobj = pymapdl_reader.read_binary(examples.fullfile)
    dofref, k, m = fobj.load_km()
    assert dofref.size
    assert k.size
    assert m.size


def test_full_sparse(sparse_full):
    str_rep = str(sparse_full)
    assert "20.1" in str_rep
    assert "MAPDL Full File" in str_rep
    assert "345" in str_rep


@skip_scipy
def test_full_sparse_k(sparse_full):
    assert isinstance(sparse_full.k, scipy.sparse.csc.csc_matrix)
    neqn = sparse_full._header["neqn"]
    assert sparse_full.k.shape == (neqn, neqn)


@skip_scipy
def test_full_sparse_m(sparse_full):
    assert isinstance(sparse_full.m, scipy.sparse.csc.csc_matrix)
    neqn = sparse_full._header["neqn"]
    assert sparse_full.m.shape == (neqn, neqn)


def test_full_sparse_dof_ref(sparse_full):
    # tests if sorted ascending
    assert (np.diff(sparse_full.dof_ref[:, 0]) >= 0).all()
    assert np.allclose(np.unique(sparse_full.dof_ref[:, 1]), [0, 1, 2])


def test_full_sparse_const(sparse_full):
    assert not sparse_full.const.any()


def test_full_load_km(sparse_full):
    dof_ref, k, m = sparse_full.load_km()
    assert not (np.diff(dof_ref[:, 0]) >= 0).all()
    neqn = sparse_full._header["neqn"]
    assert k.shape == (neqn, neqn)
    assert m.shape == (neqn, neqn)

    # make sure internal values are not overwritten
    assert (np.diff(sparse_full.dof_ref[:, 0]) >= 0).all()


def test_load_vector(sparse_full):
    assert not sparse_full.load_vector.any()


class TestPathlibFilename:
    def test_pathlib_filename_property(self, sparse_full_pathlib_full_file):
        assert isinstance(sparse_full_pathlib_full_file.pathlib_filename, pathlib.Path)

    def test_filename_property_is_string(self, sparse_full_pathlib_full_file):
        assert isinstance(sparse_full_pathlib_full_file.filename, str)

    def test_filename_setter_pathlib(self, sparse_full_pathlib_full_file):
        with pytest.raises(AttributeError):
            sparse_full_pathlib_full_file.filename = pathlib.Path("dummy2")

    def test_filename_setter_string(self, sparse_full_pathlib_full_file):
        with pytest.raises(AttributeError):
            sparse_full_pathlib_full_file.filename = "dummy2"
