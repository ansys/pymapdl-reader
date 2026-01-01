# Copyright (C) 2021 - 2026 ANSYS, Inc. and/or its affiliates.
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
