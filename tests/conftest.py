import os

import pytest
import pyvista

# Necessary for CI plotting
pyvista.OFF_SCREEN = True


@pytest.fixture(scope="session")
def mapdl(request):
    """This fixture will only be called if ``ansys.mapdl.core`` is installed."""
    from ansys.mapdl.core import launch_mapdl
    from ansys.mapdl.core.launcher import get_start_instance
    from ansys.mapdl.core.misc import get_ansys_bin

    # check if the user wants to permit pytest to start MAPDL
    # and don't allow mapdl to exit upon collection unless mapdl is local
    cleanup = get_start_instance()

    # check for a valid MAPDL install with gRPC
    valid_rver = ["211"]  # checks in this order
    EXEC_FILE = None
    for rver in valid_rver:
        if os.path.isfile(get_ansys_bin(rver)):
            EXEC_FILE = get_ansys_bin(rver)
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
