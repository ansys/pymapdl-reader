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

"""
Test loading results from plane183

Need to add ansys results for verification...

"""

import os

import numpy as np
import pytest

from ansys.mapdl import reader as pymapdl_reader

testfiles_path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def result():
    filename = os.path.join(testfiles_path, "pymapdl_182_183_42_82.rst")
    return pymapdl_reader.read_binary(filename)


def test_load(result):
    assert np.any(result.grid.cells)
    assert np.any(result.grid.points)


def test_displacement(result):
    nnum, disp = result.nodal_solution(0)
    ansys_nnum = np.load(os.path.join(testfiles_path, "prnsol_u_nnum.npy"))
    ansys_disp = np.load(os.path.join(testfiles_path, "prnsol_u.npy"))
    assert np.allclose(nnum, ansys_nnum)
    assert np.allclose(disp, ansys_disp, rtol=1e-4)  # rounding in text file


def test_stress(result):
    ansys_nnum = np.load(os.path.join(testfiles_path, "prnsol_s_nnum.npy"))
    ansys_stress = np.load(os.path.join(testfiles_path, "prnsol_s.npy"))
    nnum, stress = result.nodal_stress(0)
    mask = np.isin(nnum, ansys_nnum)
    assert np.allclose(stress[mask], ansys_stress, atol=1e-6)
