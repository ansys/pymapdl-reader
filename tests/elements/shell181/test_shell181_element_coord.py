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

Attached is the ANSYS rst file. Please consider the following expected outcomes:
element stresses in the element coordinate system:

ANSYS Version 18.2
Element Type Shell 181
Element# 118223

Top Surface
Node # Sx Sy Sz Sxy Sxz Syz
143901 2272.9 968.77 0 -222.11 0 0
143862 2291 1028.6 0 151.15 0 0
144000 3304.5 1333.3 0 171.24 0 0
143931 3286.5 1273.5 0 -202.03 0 0

Bottom Surface
Node # Sx Sy Sz Sxy Sxz Syz
143901 -2277.2 -967.94 0 243.1 0 0
143862 -2295.2 -1028.6 0 -130.16 0 0
144000 -3313 -1333.3 0 -150.24 0 0
143931 -3295 -1272.6 0 223.02 0 0

It should be noted that the average of Top and Bottom surface is what Pyansys
currently provides through result.element_stress(0) command which is in the global
coordinate system.

"""

import os

import numpy as np

from ansys.mapdl import reader as pymapdl_reader

# result for element 118223
KNOWN_RESULT_ENODE = [143901, 143862, 144000, 143931]
KNOWN_RESULT_STRESS = np.array(
    [
        [2272.9, 968.77, 0, -222.11, 0, 0],
        [2291, 1028.6, 0, 151.15, 0, 0],
        [3304.5, 1333.3, 0, 171.24, 0, 0],
        [3286.5, 1273.5, 0, -202.03, 0, 0],
        [-2277.2, -967.94, 0, 243.1, 0, 0],
        [-2295.2, -1028.6, 0, -130.16, 0, 0],
        [-3313, -1333.3, 0, -150.24, 0, 0],
        [-3295, -1272.6, 0, 223.02, 0, 0],
    ]
)

test_path = os.path.dirname(os.path.abspath(__file__))
result_file_name = os.path.join(test_path, "shell181_box.rst")


def test_shell_stress_element_cs():
    result = pymapdl_reader.read_binary(result_file_name)
    enum, stress, enode = result.element_stress(0, in_element_coord_sys=True)

    idx = np.where(enum == 118223)[0][0]
    assert np.allclose(KNOWN_RESULT_ENODE, enode[idx][:4])
    assert np.allclose(KNOWN_RESULT_STRESS, stress[idx], rtol=1e-4)
