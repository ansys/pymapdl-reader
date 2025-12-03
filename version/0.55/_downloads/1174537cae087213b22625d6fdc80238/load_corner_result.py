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
.. _ref_load_cylindrical_result:

Cylindrical Nodal Stress
~~~~~~~~~~~~~~~~~~~~~~~~
Visualize the nodal stress in the radial direction.  This is
equivalent to setting the result coordinate system to cylindrical in
MAPDL (e.g. ``RSYS, 1``).

"""

################################################################################
# Download a small result file containing the corner of a thick pipe
from ansys.mapdl.reader import examples

rst = examples.download_corner_pipe()

# obtain the cylindrical_nodal_stress
nnum, stress = rst.cylindrical_nodal_stress(0)
print(stress)

# contains results for each node in following directions
# R, THETA, Z, RTHETA, THETAZ, and RZ
print(stress.shape)

################################################################################
# plot cylindrical nodal stress in the radial direction
_ = rst.plot_cylindrical_nodal_stress(0, "R", show_edges=True, show_axes=True)

################################################################################
# plot cylindrical nodal stress in the theta direction
_ = rst.plot_cylindrical_nodal_stress(
    0, "THETA", show_edges=True, show_axes=True, add_text=False
)

################################################################################
# Plot cartesian stress in the "X" direction
_ = rst.plot_nodal_stress(0, "X", show_edges=True, show_axes=True)
