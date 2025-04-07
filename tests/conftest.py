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

import os

import pytest

from ansys.mapdl.reader.misc.checks import run_if_graphics_required

try:
    run_if_graphics_required()
    import pyvista

    # Necessary for CI plotting
    pyvista.OFF_SCREEN = True
except ImportError:
    pass

from ansys.mapdl.reader.misc.checks import are_graphics_available

skip_no_graphics = pytest.mark.skipif(
    not are_graphics_available(),
    reason="Graphic dependencies are required for this test.",
)


@pytest.fixture(scope="session")
def mapdl(request):
    """This fixture will only be called if ``ansys.mapdl.core`` is installed."""
    from ansys.mapdl.core import launch_mapdl
    from ansys.mapdl.core.launcher import get_start_instance
    from ansys.tools.path import get_mapdl_path

    # check if the user wants to permit pytest to start MAPDL
    # and don't allow mapdl to exit upon collection unless mapdl is local
    cleanup = get_start_instance()

    # check for a valid MAPDL install with gRPC
    valid_rver = ["211"]  # checks in this order
    EXEC_FILE = None
    for rver in valid_rver:
        if os.path.isfile(get_mapdl_path(rver)):
            EXEC_FILE = get_mapdl_path(rver)
            break

    return launch_mapdl(
        EXEC_FILE, override=True, cleanup_on_exit=cleanup, additional_switches="-smp"
    )


@pytest.fixture(scope="function")
def cleared(mapdl):
    mapdl.finish()
    mapdl.clear("NOSTART")  # *MUST* be NOSTART.  With START fails after 20 calls...
    mapdl.prep7()
    yield
