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

import numpy as np

from ansys.mapdl.reader.misc.checks import run_if_graphics_required

try:
    run_if_graphics_required()
    import pyvista as pv
    import vtk
except ImportError:
    pass

from conftest import skip_no_graphics

from ansys.mapdl.reader import _binary_reader

# test stress tensors from
# Sx Sy Sz Sxy Syz Sxz
stress = np.array(
    [-2.21786547, 99.05487823, -11.42874718, -4.69416809, 23.24783707, 0.4061397]
)

# known results when rotating about a vector with angle 20 degrees
# using apsg
stress_rot_x = np.array(
    [-2.21786547, 71.18732452, 16.43880463, -4.54998303, 53.31763077, -1.22385347]
)

stress_rot_y = np.array(
    [-3.03427238, 99.05487823, -10.61234027, 3.54015345, 23.45132099, -2.64919926]
)

stress_rot_z = np.array(
    [12.64614819, 84.19086457, -11.42874718, -36.1443738, 21.9847289, -7.56958209]
)


@skip_no_graphics
def test_tensor_rotation_x():
    transform = vtk.vtkTransform()
    transform.RotateX(20)
    transform.Update()
    rot_matrix = transform.GetMatrix()
    # rot_matrix.Invert()  # <-- this should not be necessary
    trans = pv.array_from_vtkmatrix(rot_matrix)

    s_test = stress.copy().reshape(1, -1)
    _binary_reader.tensor_arbitrary(s_test, trans)
    assert np.allclose(s_test, stress_rot_x)


@skip_no_graphics
def test_tensor_rotation_y():
    transform = vtk.vtkTransform()
    transform.RotateY(20)
    transform.Update()
    rot_matrix = transform.GetMatrix()
    # rot_matrix.Invert()  # <-- this should not be necessary
    trans = pv.array_from_vtkmatrix(rot_matrix)

    s_test = stress.copy().reshape(1, -1)
    _binary_reader.tensor_arbitrary(s_test, trans)
    assert np.allclose(s_test, stress_rot_y)


@skip_no_graphics
def test_tensor_rotation_z():
    transform = vtk.vtkTransform()
    transform.RotateZ(20)
    transform.Update()
    rot_matrix = transform.GetMatrix()
    # rot_matrix.Invert()  # <-- this should not be necessary
    trans = pv.array_from_vtkmatrix(rot_matrix)

    s_test = stress.copy().reshape(1, -1)
    _binary_reader.tensor_arbitrary(s_test, trans)
    assert np.allclose(s_test, stress_rot_z)
