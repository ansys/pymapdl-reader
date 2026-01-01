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

import os
import platform
import warnings

import pytest
from pyvista.plotting import system_supports_plotting

from ansys.mapdl.reader import examples

HAS_IMAGEIO = True
try:
    import imageio_ffmpeg  # noqa: F401
except ImportError:
    HAS_IMAGEIO = False

try:
    shaft = examples.download_shaft_modal()
except:
    warnings.warn("Unable to execute ``examples.download_shaft_modal``")
    shaft = None


IS_MAC = platform.system() == "Darwin"
skip_plotting = pytest.mark.skipif(
    not system_supports_plotting() or IS_MAC, reason="Requires active X Server"
)
skip_no_shaft = pytest.mark.skipif(shaft is None, reason="Requires example file")


@skip_plotting
def test_show_hex_archive():
    examples.show_hex_archive(off_screen=True)


def test_load_result():
    examples.load_result()


@skip_plotting
def test_show_displacement():
    examples.show_displacement(off_screen=True)


@skip_plotting
def test_show_stress():
    examples.show_stress(off_screen=True)


def test_load_km():
    examples.load_km()


@skip_plotting
def test_show_cell_qual():
    examples.show_cell_qual(meshtype="tet", off_screen=True)
    examples.show_cell_qual(meshtype="hex", off_screen=True)


@skip_plotting
@skip_no_shaft
@pytest.mark.skipif(not HAS_IMAGEIO, reason="Requires imageio_ffmpeg")
def test_shaft_animate(tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.mp4"))
    shaft.animate_nodal_solution(
        5,
        element_components="SHAFT_MESH",
        comp="norm",
        loop=False,
        n_frames=10,
        show_edges=True,
        off_screen=True,
        movie_filename=filename,
    )


@skip_plotting
@skip_no_shaft
def test_shaft_nodal_solution_ncomp(tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.mp4"))
    shaft.plot_nodal_solution(
        5, node_components="N_AREA_BC1", sel_type_all=False, off_screen=True
    )


@skip_plotting
@skip_no_shaft
def test_shaft_plot_screenshot(tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.png"))
    shaft.plot_nodal_solution(0, off_screen=True, screenshot=filename)
    assert os.path.isfile(filename)
