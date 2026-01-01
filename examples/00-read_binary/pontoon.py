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
.. _ref_pontoon:

Shell Static Analysis
~~~~~~~~~~~~~~~~~~~~~

Visualize a shell static analysis
"""

# download the pontoon example
from ansys.mapdl.reader import examples

pontoon = examples.download_pontoon()

###############################################################################
# Print the pontoon result
print(pontoon)

###############################################################################
# Plot the nodal displacement
pontoon.plot_nodal_solution(0, show_displacement=True, displacement_factor=100000)


###############################################################################
# print the available result types
pontoon.available_results


###############################################################################
# Plot the shell elements
pontoon.plot()

###############################################################################
# Plot the elastic strain and show exaggerated displacement
pontoon.plot_nodal_elastic_strain(
    0,
    "eqv",
    show_displacement=True,
    displacement_factor=100000,
    overlay_wireframe=True,
    lighting=False,
    add_text=False,
    show_edges=True,
)
# Note: lighting is disabled here as it's too dark
