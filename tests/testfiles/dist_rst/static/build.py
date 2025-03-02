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
These example files were built with the following script.
"""

import os

import numpy as np
import pyansys

os.environ["I_MPI_SHM_LMT"] = "shm"  # necessary on ubuntu and dmp
mapdl = pyansys.launch_mapdl(nproc=4)

mapdl.finish()
mapdl.clear()

# cylinder and mesh parameters
# torque = 100
radius = 2
h_tip = 2
height = 20
elemsize = 1.0
# pi = np.arccos(-1)
force = 100 / radius
pressure = force / (h_tip * 2 * np.pi * radius)

mapdl.prep7()
mapdl.et(1, 186)
mapdl.et(2, 154)
mapdl.r(1)
mapdl.r(2)

# Aluminum properties (or something)
mapdl.mp("ex", 1, 10e6)
mapdl.mp("nuxy", 1, 0.3)
mapdl.mp("dens", 1, 0.1 / 386.1)
mapdl.mp("dens", 2, 0)

# Simple cylinder
for i in range(4):
    mapdl.cylind(radius, "", "", height, 90 * (i - 1), 90 * i)

mapdl.nummrg("kp")

# mesh cylinder
mapdl.lsel("s", "loc", "x", 0)
mapdl.lsel("r", "loc", "y", 0)
mapdl.lsel("r", "loc", "z", 0, height - h_tip)
mapdl.lesize("all", elemsize * 2)
mapdl.mshape(0)
mapdl.mshkey(1)

mapdl.esize(elemsize)
mapdl.allsel("all")
mapdl.vsweep("ALL")
mapdl.csys(1)
mapdl.asel("s", "loc", "z", "", height - h_tip + 0.0001)
mapdl.asel("r", "loc", "x", radius)
mapdl.local(11, 1)

mapdl.csys(0)

# mesh the surface with SURF154
mapdl.aatt(2, 2, 2, 11)
mapdl.amesh("all")
mapdl.prep7()

# plot elements
# mapdl.eplot()

# Apply tangential pressure
mapdl.esel("S", "TYPE", "", 2)
mapdl.sfe("all", 2, "pres", "", pressure)

# Constrain bottom of cylinder/rod
mapdl.asel("s", "loc", "z", 0)
mapdl.nsla("s", 1)
mapdl.d("all", "all")
mapdl.allsel()

# new solution
mapdl.run("/SOLU")
mapdl.antype("static", "new")
# mapdl.eqslv('pcg', 1e-8)
mapdl.solve()
