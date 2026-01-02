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
.. _ref_load_thermal_result:

Thermal Analysis
~~~~~~~~~~~~~~~~

Visualize the result of verification manual test 33.

"""

from ansys.mapdl.reader import examples

################################################################################
# Download the result file from verification manual test case 33
vm33 = examples.download_verification_result(33)

# get nodal thermal strain for result set 1
nnum, tstrain = vm33.nodal_thermal_strain(0)

# plot nodal thermal strain for result set 11 in the X direction
vm33.plot_nodal_thermal_strain(
    10, "X", show_edges=True, lighting=True, cmap="bwr", show_axes=True
)

################################################################################
# Plot with contours

# Disable lighting and set number of colors to 10 to make an MAPDL-like plot
vm33.plot_nodal_thermal_strain(
    10,
    "X",
    show_edges=True,
    n_colors=10,
    interpolate_before_map=True,
    lighting=False,
    show_axes=True,
)
