import os
import pathlib

import numpy as np
import pytest

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader.emat import EmatFile

test_path = os.path.dirname(os.path.abspath(__file__))
testfiles_path = os.path.join(test_path, "testfiles")
emat_filename = os.path.join(testfiles_path, "file.emat")


@pytest.fixture(scope="module")
def emat():
    emat_bin = pymapdl_reader.read_binary(emat_filename)
    assert isinstance(emat_bin, EmatFile)
    return emat_bin


@pytest.fixture(scope="module")
def emat_pathlib():
    emat_bin = EmatFile(pathlib.Path(emat_filename))
    return emat_bin


def test_load_element(emat):
    dof_idx, element_data = emat.read_element(0)
    assert "stress" in element_data
    assert "mass" in element_data


def test_global_applied_force(emat):
    force = emat.global_applied_force()
    assert np.allclose(force, 0)


def test_eeqv(emat):
    assert np.allclose(np.sort(emat.eeqv), emat.enum)


def test_neqv(emat):
    assert np.allclose(np.sort(emat.neqv), emat.nnum)


class TestPathlibFilename:
    def test_pathlib_filename_property(self, emat_pathlib):
        assert isinstance(emat_pathlib.pathlib_filename, pathlib.Path)

    def test_filename_property_is_string(self, emat_pathlib):
        assert isinstance(emat_pathlib.filename, str)

    def test_filename_setter_pathlib(self, emat_pathlib):
        with pytest.raises(AttributeError):
            emat_pathlib.filename = pathlib.Path("dummy2")

    def test_filename_setter_string(self, emat_pathlib):
        with pytest.raises(AttributeError):
            emat_pathlib.filename = "dummy2"
