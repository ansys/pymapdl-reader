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

"""
Test loading results from shell181 elements

Reading the stress values from these elements requires more overhead as
the component stresses written to the binary file are relative to the element's
coordinates and not the global coordinates.


Element solution results from ANSYS
###############################################################################
 PRINT S    ELEMENT SOLUTION PER ELEMENT

 ***** POST1 ELEMENT NODAL STRESS LISTING *****

  LOAD STEP=     1  SUBSTEP=     1
   TIME=    1.0000      LOAD CASE=   0
  SHELL RESULTS FOR TOP/BOTTOM ALSO MID WHERE APPROPRIATE

  THE FOLLOWING X,Y,Z VALUES ARE IN GLOBAL COORDINATES


  ELEMENT=       1        SHELL181
    NODE    SX          SY          SZ          SXY         SYZ         SXZ
       2  0.17662E-07  79.410     -11.979    -0.11843E-02  4.8423    -0.72216E-04
       1  0.20287E-07  91.212      27.364    -0.13603E-02  4.8423    -0.72216E-04
       4  0.20287E-07  91.212      27.364    -0.13603E-02 -4.8423     0.72216E-04
       3  0.17662E-07  79.410     -11.979    -0.11843E-02 -4.8423     0.72216E-04
       2 -0.17662E-07 -79.410      11.979     0.11843E-02 -4.8423     0.72216E-04
       1 -0.20287E-07 -91.213     -27.364     0.13603E-02 -4.8423     0.72216E-04
       4 -0.20287E-07 -91.213     -27.364     0.13603E-02  4.8423    -0.72216E-04
       3 -0.17662E-07 -79.410      11.979     0.11843E-02  4.8423    -0.72216E-04
###############################################################################


Nodal Solution results from ANSYS
###############################################################################
 PRINT S    NODAL SOLUTION PER NODE

  ***** POST1 NODAL STRESS LISTING *****
  PowerGraphics Is Currently Enabled

  LOAD STEP=     1  SUBSTEP=     1
   TIME=    1.0000      LOAD CASE=   0
  SHELL NODAL RESULTS ARE AT TOP/BOTTOM FOR MATERIAL   1

  THE FOLLOWING X,Y,Z VALUES ARE IN GLOBAL COORDINATES

    NODE    SX          SY          SZ          SXY         SYZ         SXZ
       1  0.20287E-07  91.212      27.364    -0.13603E-02  4.8423    -0.72216E-04
       1 -0.20287E-07 -91.213     -27.364     0.13603E-02 -4.8423     0.72216E-04
       2  0.17662E-07  79.410     -11.979    -0.11843E-02  4.8423    -0.72216E-04
       2 -0.17662E-07 -79.410      11.979     0.11843E-02 -4.8423     0.72216E-04
       3  0.17662E-07  79.410     -11.979    -0.11843E-02 -4.8423     0.72216E-04
       3 -0.17662E-07 -79.410      11.979     0.11843E-02  4.8423    -0.72216E-04
       4  0.20287E-07  91.212      27.364    -0.13603E-02 -4.8423     0.72216E-04
       4 -0.20287E-07 -91.213     -27.364     0.13603E-02  4.8423    -0.72216E-04

  ***** POST1 NODAL STRESS LISTING *****

  LOAD STEP=     1  SUBSTEP=     1
   TIME=    1.0000      LOAD CASE=   0
  SHELL NODAL RESULTS ARE AT TOP/BOTTOM FOR MATERIAL   4

  THE FOLLOWING X,Y,Z VALUES ARE IN GLOBAL COORDINATES

    NODE    SX          SY          SZ          SXY         SYZ         SXZ

 MINIMUM VALUES
 NODE          1           1           1           1           1           1
 VALUE  -0.20287E-07 -91.213     -27.364    -0.13603E-02 -4.8423    -0.72216E-04

 MAXIMUM VALUES
 NODE          1           1           1           1           1           1
 VALUE   0.20287E-07  91.212      27.364     0.13603E-02  4.8423     0.72216E-04
###############################################################################


"""

import os

import numpy as np
import pytest

from ansys.mapdl import reader as pymapdl_reader

ANSYS_ELEM = [
    [0.17662e-07, 79.410, -11.979, -0.11843e-02, 4.8423, -0.72216e-04],
    [0.20287e-07, 91.212, 27.364, -0.13603e-02, 4.8423, -0.72216e-04],
    [0.20287e-07, 91.212, 27.364, -0.13603e-02, -4.8423, 0.72216e-04],
    [0.17662e-07, 79.410, -11.979, -0.11843e-02, -4.8423, 0.72216e-04],
]

ANSYS_NODE = [
    [0.20287e-07, 91.212, 27.364, -0.13603e-02, 4.8423, -0.72216e-04],
    [0.17662e-07, 79.410, -11.979, -0.11843e-02, 4.8423, -0.72216e-04],
    [0.17662e-07, 79.410, -11.979, -0.11843e-02, -4.8423, 0.72216e-04],
    [0.20287e-07, 91.212, 27.364, -0.13603e-02, -4.8423, 0.72216e-04],
]


@pytest.fixture(scope="module")
def result():
    test_path = os.path.dirname(os.path.abspath(__file__))
    return pymapdl_reader.read_binary(os.path.join(test_path, "shell181.rst"))


def test_load(result):
    assert np.any(result.grid.cells)
    assert np.any(result.grid.points)


def test_element_stress(result):
    _, element_stress, _ = result.element_stress(0)
    element0 = element_stress[0]

    # ansys prints both positive and negative component values
    if np.sign(element0[0][0]) != np.sign(ANSYS_ELEM[0][0]):
        element0 *= -1

    # wide atol limits considering the 5 sigfig from ASCII tables
    assert np.allclose(element0, np.array(ANSYS_ELEM), atol=1e-6)


def test_nodal_stress(result):
    _, stress = result.nodal_stress(0)
    if np.sign(stress[0][0]) != np.sign(ANSYS_NODE[0][0]):
        stress *= -1

    # wide atol limits considering the 5 sigfig from ASCII tables
    assert np.allclose(stress, np.array(ANSYS_NODE), atol=1e-6)
