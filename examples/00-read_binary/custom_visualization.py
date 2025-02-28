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
.. _ref_custom_visualization:

Custom Scalar Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Display custom scalars using an existing mesh.

"""

import numpy as np

from ansys.mapdl.reader import examples

# Download an example shaft modal analysis result file
shaft = examples.download_shaft_modal()

###############################################################################
# Each result file contains both a ``mesh`` property and a ``grid``
# property.  The ``mesh`` property can be through as the MAPDL
# representation of the FEM while the ``grid`` property can be through
# of the Python visualizing property used to plot within Python.

print("shaft.mesh:\n", shaft.mesh)
print("-" * 79)
print("shaft.grid:\n", shaft.grid)

###############################################################################
# Plotting
# ~~~~~~~~
#
# The grid instance is a :class:`pyvista.UnstructuredGrid` part of the `pyvista
# <https://docs.pyvista.org/>`_ library.  This class allows for advanced
# plotting using VTK in just a few lines of code.  For example, you can plot
# the underlying mesh with:

shaft.grid.plot(color="w", smooth_shading=True)

###############################################################################
# Plotting Node Scalars
# ~~~~~~~~~~~~~~~~~~~~~
#
# If you point-wise or cell-wise scalars (nodes and elements in FEA),
# you can plot these scalars by setting the ``scalars=`` parameter.
# Here, I'm simply using the x location of the nodes to color the
# mesh.
#
# It follows that you can use any set of scalars provided that it
# matches the number of nodes in the unstructured grid or the number
# of cells in the unstructured grid.  Here, we're plotting node values.

x_scalars = shaft.grid.points[:, 0]
shaft.grid.plot(scalars=x_scalars, smooth_shading=True)


###############################################################################
# Plotting With Missing Values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# If you do not have values for every node (for example, the midside
# nodes), you can leave these values as NAN and the plotter will take
# care of plotting only the real values.
#
# For example, if you have calculated strain scalars that are only
# available at certain nodes, you can still plot those.  This example
# just nulls out the first 2000 nodes to be able to visualize the
# missing values.  If your are just missing midside values, your plot
# will not show the missing values since `pyvista` only plots the edge
# nodes.

pontoon = examples.download_pontoon()
nnum, strain = pontoon.nodal_elastic_strain(0)

scalars = strain[:, 0]
scalars[:2000] = np.nan  # here, we simulate unknown values

pontoon.grid.plot(scalars=scalars, show_edges=True, lighting=False)
