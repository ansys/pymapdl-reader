"""PyMAPDL-Reader example files"""
import inspect
import os
import sys

import numpy as np
import pyvista as pv

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples

# get location of this folder and the example files
dir_path = os.path.dirname(os.path.realpath(__file__))
rstfile = os.path.join(dir_path, "file.rst")
hexarchivefile = os.path.join(dir_path, "HexBeam.cdb")
tetarchivefile = os.path.join(dir_path, "TetBeam.cdb")
fullfile = os.path.join(dir_path, "file.full")
sector_archive_file = os.path.join(dir_path, "sector.cdb")


def run_all(run_ansys=False):
    """
    Runs all the functions within this module except for the ansys
    tests.

    """
    testfunctions = []
    for name, obj in inspect.getmembers(sys.modules[__name__]):
        if inspect.isfunction(obj) and name != "run_all":
            if "ansys" in name and not run_ansys:
                continue
            testfunctions.append(obj)

    # run all the functions
    any(f() for f in testfunctions)


def show_hex_archive(off_screen=None):
    """Display a hex beam mesh"""
    # Load an archive file
    archive = pymapdl_reader.Archive(hexarchivefile)
    archive.plot(off_screen=off_screen, color="w", show_edges=True)
    assert archive.grid.n_points
    assert archive.grid.n_cells


def load_result():
    """Loads a result file and prints out the displacement of all the nodes from
    a modal analysis.
    """

    # Load result file
    result = pymapdl_reader.read_binary(rstfile)
    assert result.nsets == 6
    assert len(result.mesh.nnum) == 321
    print("Loaded result file with {:d} result sets".format(result.nsets))
    print("Contains {:d} nodes".format(len(result.mesh.nnum)))

    # display result
    nnum, disp = result.nodal_solution(0)

    print("Nodal displacement for nodes 30 to 40 is:")

    for i in range(29, 40):
        node = result.mesh.nnum[i]
        x = disp[i, 0]
        y = disp[i, 1]
        z = disp[i, 2]
        print("{:2d}  {:10.6f}   {:10.6f}   {:10.6f}".format(node, x, y, z))


def show_displacement(off_screen=None):
    """Load and plot 1st bend of a hexahedral beam"""

    # get location of this file
    fobj = pymapdl_reader.read_binary(rstfile)

    print("Displaying ANSYS Mode 1")
    fobj.plot_nodal_solution(
        0, label="Displacement", off_screen=off_screen, n_colors=9, show_edges=True
    )


def show_stress(off_screen=None):
    """Load and plot 1st bend of a hexahedral beam"""

    # get location of this file
    result = pymapdl_reader.read_binary(rstfile)

    print("Displaying node averaged stress in x direction for Mode 6")
    result.plot_nodal_stress(5, "x", off_screen=off_screen, n_colors=9)


def load_km():
    """Loads m and k matrices from a full file"""

    # Create file reader object
    fobj = pymapdl_reader.read_binary(fullfile)
    dofref, k, m = fobj.load_km()

    # print results
    ndim = k.shape[0]
    print("Loaded {:d} x {:d} mass and stiffness matrices".format(ndim, ndim))
    print("\t k has {:d} entries".format(k.indices.size))
    print("\t m has {:d} entries".format(m.indices.size))

    # compute natural frequencies if installed
    try:
        from scipy import sparse
        from scipy.sparse import linalg
    except ImportError:
        return

    k += sparse.triu(k, 1).T
    m += sparse.triu(m, 1).T
    k += sparse.diags(np.random.random(k.shape[0]) / 1e20, shape=k.shape)

    # Solve
    w, v = linalg.eigsh(k, k=20, M=m, sigma=10000)

    # System natural frequencies
    freq = np.real(w) ** 0.5 / (2 * np.pi)

    print("First four natural frequencies:")
    for i in range(4):
        print("{:.3f} Hz".format(freq[i]))

    known_result = np.array(
        [
            1283.20036921,
            1283.20036921,
            5781.97486169,
            6919.39887714,
            6919.39887714,
            10172.61497694,
            16497.85701889,
            16497.85701889,
            17343.9939669,
            27457.18472747,
            27457.18472747,
            28908.52552073,
            30326.16886062,
            39175.76412419,
            39175.76412419,
            40503.70406456,
            49819.91597612,
            51043.03965541,
            51043.03965541,
            52193.86143879,
        ]
    )
    assert np.allclose(freq, known_result)


def solve_km():
    """Load and solves a mass and stiffness matrix from an ansys full file"""
    try:
        from scipy import sparse
        from scipy.sparse import linalg
    except ImportError:
        print("scipy not installed, aborting")
        return

    # load the mass and stiffness matrices
    full = pymapdl_reader.read_binary(examples.fullfile)
    dofref, k, m = full.load_km(sort=True)

    # make symmetric
    k += sparse.triu(k, 1).T
    m += sparse.triu(m, 1).T
    k += sparse.diags(np.random.random(k.shape[0]) / 1e20, shape=k.shape)

    # Solve
    w, v = linalg.eigsh(k, k=20, M=m, sigma=1000)

    # System natural frequencies
    f = (np.real(w)) ** 0.5 / (2 * np.pi)

    # %% Plot result

    # Get the 4th mode shape
    full_mode_shape = v[:, 3]  # x, y, z displacement for each node

    # reshape and compute the normalized displacement
    disp = full_mode_shape.reshape((-1, 3))
    n = (disp * disp).sum(1) ** 0.5
    n /= n.max()  # normalize

    # load an archive file and create a vtk unstructured grid
    archive = pymapdl_reader.Archive(examples.hexarchivefile)
    grid = archive.grid

    # Fancy plot the displacement
    pl = pv.Plotter()

    # add two meshes to the plotting class
    pl.add_mesh(grid.copy(), color="w", style="wireframe")
    pl.add_mesh(
        grid,
        scalars=n,
        stitle="Normalized\nDisplacement",
        flip_scalars=True,
        cmap="viridis",
    )
    # Update the coordinates by adding the mode shape to the grid
    pl.update_coordinates(grid.points + disp / 80, render=False)
    pl.add_text("Cantliver Beam 4th\nMode Shape at\n{:.4f}".format(f[3]), font_size=30)
    pl.plot()


def show_cell_qual(meshtype="tet", **kwargs):
    """Displays minimum scaled jacobian of a sample mesh.

    For an ANSYS analysis to run properly, the minimum scaled Jacobian
    of each cell must be greater than 0.

    Parameters
    ----------
    meshtype : str, optional
        Set to 'hex' to display cell quality of a hexahedral meshed beam.

    """

    # load archive file and parse for subsequent FEM queries
    if meshtype == "hex":
        archive = pymapdl_reader.Archive(hexarchivefile)
    else:
        archive = pymapdl_reader.Archive(tetarchivefile)

    # get cell quality
    qual = archive.quality
    assert np.all(qual > 0)

    # plot cell quality
    archive.grid.plot(
        scalars=qual,
        scalar_bar_args={"title": "Cell Minimum Scaled\nJacobian"},
        cmap="bwr",
        show_edges=True,
        rng=[0, 1],
        **kwargs,
    )
