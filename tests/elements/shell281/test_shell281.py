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

import pytest

from ansys.mapdl import reader as pymapdl_reader

path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def rst():
    rst_file = os.path.join(path, "file.rst")
    return pymapdl_reader.read_binary(rst_file)


def test_materials(rst):
    materials = rst.materials
    material = materials[1]
    material["EX"] = 40000000000
    material["EY"] = 10000000000
    material["EZ"] = 10000000000
    material["PRXY"] = 0.3
    material["PRYZ"] = 0.3
    material["PRXZ"] = 0.3
    material["GXY"] = 5000000000
    material["GYZ"] = 5000000000
    material["GXZ"] = 5000000000


def test_sections(rst):
    sections = rst.section_data
    assert isinstance(sections, dict)
    # assert isinstance(sections[3], np.ndarray)

    # TODO: add known result types for the section data
