import os
import pathlib

import numpy as np
import pytest
import scipy

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples
from ansys.mapdl.reader.full import FullFile

test_path = os.path.dirname(os.path.abspath(__file__))
testfiles_path = os.path.join(test_path, "testfiles")


@pytest.fixture()
def sparse_full_pathlib_full_file():
    filename = os.path.join(testfiles_path, "sparse.full")
    return FullFile(pathlib.Path(filename))


@pytest.fixture()
def sparse_full():
    filename = os.path.join(testfiles_path, "sparse.full")
    return pymapdl_reader.read_binary(filename)


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


def test_full_sparse_k(sparse_full):
    assert isinstance(sparse_full.k, scipy.sparse.csc.csc_matrix)
    neqn = sparse_full._header["neqn"]
    assert sparse_full.k.shape == (neqn, neqn)


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
