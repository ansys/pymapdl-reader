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

import numpy as np
import pytest

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples

test_path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def result():
    return pymapdl_reader.read_binary(examples.rstfile)


@pytest.fixture(scope="module")
def archive():
    return pymapdl_reader.Archive(examples.hexarchivefile)


def test_geometry_elements(result, archive):
    r_elem = np.array(result.mesh.elem)[result._sidx_elem]
    assert np.allclose(r_elem, archive.elem)


def test_geometry_nodes(result, archive):
    assert np.allclose(result.mesh.nodes[:, :3], archive.nodes)


def test_geometry_nodenum(result, archive):
    assert np.allclose(result.mesh.nnum, archive.nnum)


def test_results_displacement(result):
    textfile = os.path.join(test_path, "prnsol_u.txt")
    nnum, r_values = result.nodal_solution(0)
    a_values = np.loadtxt(textfile, skiprows=2)[:, 1:4]
    assert np.allclose(r_values, a_values)


def test_results_stress(result):
    _, r_values = result.nodal_stress(0)
    textfile = os.path.join(test_path, "prnsol_s.txt")
    a_values = np.loadtxt(textfile, skiprows=2)[:, 1:]

    # ignore nan
    nanmask = ~np.isnan(r_values).any(1)
    assert np.allclose(r_values[nanmask], a_values, atol=1e-1)


def test_results_pstress(result):
    r_nnum, r_values = result.principal_nodal_stress(0)
    textfile = os.path.join(test_path, "prnsol_s_prin.txt")
    a_values = np.loadtxt(textfile, skiprows=2)[:, 1:]

    # ignore nan
    nanmask = ~np.isnan(r_values).any(1)
    assert np.allclose(r_values[nanmask], a_values, atol=100)
