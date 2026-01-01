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
.. _ref_sector_model:

Cyclic Model Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~
Visualize and animate a full cyclic model.  This model is based on the
jetcat rotor.

First, load the rotor.  Notice how printing the rotor class reveals
the details of the rotor result file.
"""

# sphinx_gallery_thumbnail_number = 2
from ansys.mapdl.reader import examples

rotor = examples.download_sector_modal()
print(rotor)

###############################################################################
# Plot the rotor and rotor sectors
#
# Note that additional keyword arguments can be passed to the plotting
# functions of ``pymapdl-reader``.  See ``help(pyvista.plot`` for the
# documentation on all the keyword arguments.
rotor.plot_sectors(cpos="xy", smooth_shading=True)
rotor.plot()


###############################################################################
# Plot nodal displacement for result 21.
#
# Note that pymapdl-reader uses 0 based cumulative indexing.  You could also
# use the (load step, sub step) ``(4, 3)``.
rotor.plot_nodal_displacement(
    20, show_displacement=True, displacement_factor=0.001, overlay_wireframe=True
)  # same as (2, 4)


###############################################################################
# Animate Mode 21
# ~~~~~~~~~~~~~~~
# Disable movie_filename and increase n_frames for a smoother plot
rotor.animate_nodal_solution(
    20,
    loop=False,
    movie_filename="rotor_mode.gif",
    background="w",
    displacement_factor=0.001,
    add_text=False,
    n_frames=30,
)
